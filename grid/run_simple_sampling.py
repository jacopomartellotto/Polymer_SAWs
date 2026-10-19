#!/usr/bin/env python3
"""Simple-sampling rejection method for 2D self-avoiding walks.

Generate non-backtracking random walks on the square lattice and reject the
whole trial at the first repeated site. The first step chooses four directions;
subsequent steps choose the three directions other than the immediate reverse.
Each N-step proposal has probability 1 / (4 * 3**(N-1)), so conditioning on
survival still samples N-step SAWs uniformly. All other revisits cause rejection.
"""

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd

# N conta i passi: un cammino completo contiene N+1 siti, origine inclusa.
# La spaziatura del reticolo vale 1 e non ci sono energie di interazione.
DEFAULT_LENGTHS = "4,5,6,7,8,9,10,12,14,16,18,20,22,25,28,30"
PROPOSAL_SCHEME = "non_backtracking"

# Tentativi aggiuntivi: aggiungo statistica dove i sopravvissuti sono più rari.
THESIS_BUDGETS = {20: 2_000_000, 22: 4_000_000, 25: 12_000_000,
                  28: 30_000_000, 30: 60_000_000}


#comandi CLI 
def parse_lengths(text):
    # Converte la lista di lunghezze da riga di comando in interi.
    values = []
    for part in text.split(","):
        part = part.strip()
        if part:
            values.append(int(part))
    if not values:
        raise ValueError("at least one length is required")
    if any(length < 1 for length in values):
        raise ValueError("walk lengths must be positive")
    return values


#Genero i cammini
def generate_simple_saw(length, rng):
    # Escludere il ritorno immediato evita rigetti certi, senza favorire
    # alcuni SAW rispetto ad altri: tutte le proposte complete sono equiprobabili.
    x = 0
    y = 0
    # Memorizziamo soltanto i siti occupati
    visited = {(0, 0)}
    directions = ((1, 0), (-1, 0), (0, 1), (0, -1))
    previous_direction = None

    for step in range(1, length + 1):
        if previous_direction is None:
            # Primo passo: quattro direzioni, ciascuna con probabilita' 1/4.
            direction = int(rng.integers(0, 4))
        else:
            # Le coppie di indici 0/1 e 2/3 sono opposte: XOR 1 trova
            # il ritorno immediato. Si salta SOLO questo indice, non gli
            # altri siti occupati, che devono ancora causare il rigetto.
            forbidden = previous_direction ^ 1
            direction = int(rng.integers(0, 3))
            direction += direction >= forbidden
        dx, dy = directions[direction]
        x += dx
        y += dy
        pos = (x, y)
        if pos in visited:
            # Appena torna su un sito gia' visitato, la camminata non e' self-avoiding.
            # Il passo che causa la collisione non entra in L: sono validi soltanto i precedenti step-1 passi.
            return None, step - 1, True
        visited.add(pos)
        previous_direction = direction

    # L'origine resta (0,0): Re2 e' il quadrato della norma del punto finale.
    # Un successo contribuisce con L=N, quindi con frazione normalizzata 1.
    return x * x + y * y, length, False


def simulate_length(length, walks, rng):
    # Ripete il simple sampling per una lunghezza e riassume le statistiche.
    r2_values = []
    free_lengths = np.empty(walks, dtype=np.int32)
    failures = 0

    for i in range(walks):
        # Re2 si misura SOLO sui cammini completi. Tutti i tentativi,
        # riusciti e falliti, contribuiscono invece alle statistiche di L
        # e alle probabilita' di successo/fallimento.
        r2, free_length, failed = generate_simple_saw(length, rng)
        free_lengths[i] = free_length
        if failed:
            failures += 1
        else:
            r2_values.append(float(r2))

    survivors = walks - failures
    survival = survivors / walks
    failure = failures / walks
    mean_free_path = float(np.mean(free_lengths))
    free_path_sd = float(np.std(free_lengths, ddof=1)) if walks > 1 else math.nan
    free_path_sem = free_path_sd / math.sqrt(walks) if walks > 1 else math.nan
    # N e' costante in questa run: media(L/N) = media(L)/N. Nel Pivot
    # la coda ha lunghezza variabile e questa semplificazione non vale.
    mean_free_path_fraction = mean_free_path / length if length > 0 else math.nan

    if survivors > 1:
        # Deviazione standard e SEM hanno senso solo con almeno due sopravvissuti.
        r2_array = np.asarray(r2_values, dtype=float)
        mean_r2 = float(np.mean(r2_array))
        sd_r2 = float(np.std(r2_array, ddof=1))
        # SD: dispersione delle singole configurazioni. SEM: errore della
        # media su tentativi indipendenti; usa i sopravvissuti, non i tentativi.
        sem_r2 = sd_r2 / math.sqrt(survivors)
    elif survivors == 1:
        mean_r2 = r2_values[0]
        sd_r2 = math.nan
        sem_r2 = math.nan
    else:
        mean_r2 = math.nan
        sd_r2 = math.nan
        sem_r2 = math.nan

    # Errore standard binomiale stimato da tentativi indipendenti.
    # Successo e fallimento sono complementari e hanno lo stesso SEM;
    # questa formula non e' un intervallo di confidenza esatto ai bordi 0 e 1.
    survival_sem = math.sqrt(max(survival * (1.0 - survival), 0.0) / walks)

    return {
        "proposal_scheme": PROPOSAL_SCHEME,
        "N": length,
        "walks": walks,
        "survivors": survivors,
        "failures": failures,
        "survival": survival,
        "survival_sem": survival_sem,
        "failure": failure,
        "failure_sem": survival_sem,
        "Re2_mean": mean_r2,
        "Re2_sd": sd_r2,
        "Re2_sem": sem_r2,
        "mean_free_path": mean_free_path,
        "mean_free_path_sd": free_path_sd,
        "mean_free_path_sem": free_path_sem,
        "mean_free_path_fraction": mean_free_path_fraction,
        "mean_free_path_fraction_sem": free_path_sem / length if length > 0 else math.nan,
    }


def vectorized_moments(length, count, rng, batch_size=100_000):
    """Return moments of independent non-backtracking rejection trials."""
    # Stesso campionamento del metodo scalare, eseguito per blocchi.
    # Conserva: tentativi, sopravvissuti, somma(Re2), somma(Re2^2),
    # somma(L), somma(L^2). Bastano per riunire medie e varianze.
    totals = np.zeros(6, dtype=np.int64)
    # Codifica un sito come un intero su un quadrato di lato 2N+1.
    # In N passi dall'origine le coordinate restano fra -N e N:
    # la codifica e' univoca e non introduce bordi fisici nella simulazione.
    stride = 2 * length + 1
    origin = length * stride + length
    deltas = np.array([stride, -stride, 1, -1], dtype=np.int32)
    for start in range(0, count, batch_size):
        size = min(batch_size, count - start)
        # Ogni riga e' un tentativo ancora vivo; le colonne sono i siti visitati.
        history = np.full((size, 1), origin, dtype=np.int32)
        previous_directions = None
        length_sum = 0
        length_sq_sum = 0
        for step in range(1, length + 1):
            if previous_directions is None:
                choices = rng.integers(0, 4, size=len(history))
            else:
                forbidden = previous_directions ^ 1
                choices = rng.integers(0, 3, size=len(history))
                choices += choices >= forbidden
            candidate = history[:, -1] + deltas[choices]
            # Confronto con tutta la storia della stessa camminata.
            hit = (history == candidate[:, None]).any(axis=1)
            failures = int(hit.sum())
            length_sum += failures * (step - 1)
            length_sq_sum += failures * (step - 1) ** 2
            # I falliti vengono eliminati dal blocco dopo aver contato L;
            # i superstiti proseguono con nuove direzioni uniformi.
            history = np.column_stack([history[~hit], candidate[~hit]])
            # Conserva la direzione precedente dei soli tentativi sopravvissuti.
            previous_directions = choices[~hit]
            if len(history) == 0:
                break
        survivors = len(history)
        length_sum += survivors * length
        length_sq_sum += survivors * length * length
        endpoint = history[:, -1]
        r2 = ((endpoint // stride - length) ** 2
              + (endpoint % stride - length) ** 2).astype(np.int64)
        totals += [size, survivors, r2.sum(), (r2 * r2).sum(),
                   length_sum, length_sq_sum]
    return totals


def refine_thesis_data(base_csv, force=False, out_dir=None, budgets=None,
                       refinement_seed=2026092000):
    """Apply the fixed high-statistics thesis protocol and pool all moments."""
    root = Path(__file__).resolve().parents[1]
    out = (Path(out_dir) if out_dir else
           root / "outputs/data/simple_sampling/refined")
    out.mkdir(parents=True, exist_ok=True)
    # Si parte sempre dai dati di base, non dal precedente file raffinato:
    # rieseguire il pooling non deve contare due volte le stesse osservazioni.
    original = pd.read_csv(base_csv).set_index("N")
    if ("proposal_scheme" not in original
            or not original.proposal_scheme.eq(PROPOSAL_SCHEME).all()):
        raise RuntimeError("Base data use a different proposal scheme; regenerate them")
    combined = original.copy()
    budgets = THESIS_BUDGETS if budgets is None else budgets
    protocol = {
        "budgets": {str(n): count for n, count in budgets.items()},
        "base_seed": refinement_seed,
        "batch_size": 100_000,
        "proposal_scheme": PROPOSAL_SCHEME,
        "method": "independent non-backtracking walks, direct rejection at first intersection",
    }
    (out / "manifest.json").write_text(json.dumps(protocol, indent=2))

    for length, budget in budgets.items():
        moment_path = out / f"moments_N{length}.json"
        # I momenti salvati permettono di riprendere il protocollo senza
        # rigenerare i blocchi gia' completati. force ne richiede il ricalcolo.
        if moment_path.exists() and not force:
            saved = json.loads(moment_path.read_text())
            if (saved.get("proposal_scheme") != PROPOSAL_SCHEME
                    or saved.get("seed") != refinement_seed + length):
                raise RuntimeError(
                    f"Saved moments for N={length} use a different scheme or seed; "
                    "use --force-refinement to regenerate them"
                )
            moments = np.asarray(saved["moments"], dtype=np.int64)
            if int(moments[0]) != budget:
                raise RuntimeError(f"Saved moments for N={length} use a different budget")
        else:
            started = time.perf_counter()
            rng = np.random.default_rng(refinement_seed + length)
            moments = vectorized_moments(length, budget, rng, protocol["batch_size"])
            moment_path.write_text(json.dumps({
                "proposal_scheme": PROPOSAL_SCHEME,
                "moments": moments.tolist(),
                "seed": refinement_seed + length,
                "seconds": time.perf_counter() - started,
            }, indent=2))

        old = original.loc[length]
        trials = int(old.walks) + int(moments[0])
        survivors = int(old.survivors) + int(moments[1])
        # Pooling pesato con le numerosita' effettive: non si fa la media
        # semplice delle medie. Ricostruiamo la somma dei quadrati dalla SD
        # campionaria usando sum(x^2) = (n-1)*SD^2 + n*media^2.
        r_sum = old.Re2_mean * old.survivors + moments[2]
        r_sq_sum = ((old.survivors - 1) * old.Re2_sd ** 2
                    + old.survivors * old.Re2_mean ** 2 + moments[3])
        r_mean = r_sum / survivors
        r_sd = math.sqrt((r_sq_sum - survivors * r_mean ** 2) / (survivors - 1))
        # Per L la numerosita' e' il totale dei tentativi; per Re2 era
        # il numero dei sopravvissuti. Confonderli falserebbe i SEM.
        free_sum = old.mean_free_path * old.walks + moments[4]
        free_sq_sum = ((old.walks - 1) * old.mean_free_path_sd ** 2
                       + old.walks * old.mean_free_path ** 2 + moments[5])
        free_mean = free_sum / trials
        free_sd = math.sqrt((free_sq_sum - trials * free_mean ** 2) / (trials - 1))
        survival = survivors / trials
        binomial_sem = math.sqrt(survival * (1 - survival) / trials)
        values = {
            "walks": trials, "survivors": survivors, "failures": trials - survivors,
            "survival": survival, "survival_sem": binomial_sem,
            "failure": 1 - survival, "failure_sem": binomial_sem,
            "Re2_mean": r_mean, "Re2_sd": r_sd, "Re2_sem": r_sd / math.sqrt(survivors),
            "mean_free_path": free_mean, "mean_free_path_sd": free_sd,
            "mean_free_path_sem": free_sd / math.sqrt(trials),
            "mean_free_path_fraction": free_mean / length,
            "mean_free_path_fraction_sem": free_sd / math.sqrt(trials) / length,
            "reliable": int(survivors >= 100), "seconds": math.nan,
        }
        for column, value in values.items():
            combined.loc[length, column] = value
        print(f"N={length}: pooled {survivors} survivors, Re2={r_mean:.5f}")

    combined.to_csv(out / "grid_refined.csv")
    print(f"Wrote {out / 'grid_refined.csv'}")


def main():
    # Configura il batch di lunghezze, il numero di cammini ordinari provati
    # e il file CSV in cui salvare il baseline del simple sampling.
    parser = argparse.ArgumentParser(description="Simple-sampling rejection method for grid SAWs.")
    parser.add_argument("--lengths", default=DEFAULT_LENGTHS, help="comma-separated N values")
    parser.add_argument("--walks", type=int, default=1_000_000, help="non-backtracking trials per length")
    parser.add_argument("--seed", type=int, default=20260619)
    parser.add_argument("--min-survivors", type=int, default=100)
    parser.add_argument("--out", default=os.path.join(
        "outputs", "data", "simple_sampling", "base.csv"))
    parser.add_argument("--refinement-out", default=os.path.join(
        "outputs", "data", "simple_sampling", "refined"))
    parser.add_argument("--protocol", help="central simulation_protocol.json")
    parser.add_argument("--thesis", action="store_true",
                        help="also run the fixed high-statistics thesis refinement")
    parser.add_argument("--force-refinement", action="store_true",
                        help="recompute saved high-statistics sufficient statistics")
    args = parser.parse_args()

    refinement_budgets = THESIS_BUDGETS
    refinement_seed = 2026092000
    expected_lengths = parse_lengths(DEFAULT_LENGTHS)
    expected_walks = 1_000_000
    expected_seed = 20260619
    if args.protocol:
        configured = json.loads(Path(args.protocol).read_text(encoding="utf-8"))["simple_sampling"]
        if configured.get("proposal_scheme") != PROPOSAL_SCHEME:
            raise SystemExit("Configured proposal scheme does not match this sampler")
        expected_lengths = [int(n) for n in configured["lengths"]]
        expected_walks = int(configured["base_walks_per_length"])
        expected_seed = int(configured["base_seed"])
        refinement_seed = int(configured["refinement_seed"])
        refinement_budgets = {
            int(n): int(count) for n, count in configured["refinement_walks"].items()
        }

    lengths = parse_lengths(args.lengths)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    # Un unico flusso casuale attraversa le lunghezze in ordine. Stesso seed,
    # protocollo e ambiente consentono di riprodurre il flusso; cambiare
    # ordine, budget o implementazione puo' cambiare i singoli dati.
    rng = np.random.default_rng(args.seed)

    rows = []
    start_all = time.perf_counter()
    for length in lengths:
        # Ogni riga del CSV finale corrisponde a una lunghezza N.
        start = time.perf_counter()
        row = simulate_length(length, args.walks, rng)
        row["seconds"] = time.perf_counter() - start
        # Flag operativo basato su una soglia di sopravvissuti, non una
        # dimostrazione della precisione o della validita' di un fit.
        row["reliable"] = int(row["survivors"] >= args.min_survivors and math.isfinite(row["Re2_mean"]))
        rows.append(row)
        print(
            f"N={length:4d} survivors={row['survivors']:8d}/{args.walks} "
            f"P_surv={row['survival']:.6g} Re2={row['Re2_mean']:.6g} "
            f"reliable={row['reliable']}"
        )

    fieldnames = [
        "proposal_scheme",
        "N",
        "walks",
        "survivors",
        "failures",
        "survival",
        "survival_sem",
        "failure",
        "failure_sem",
        "Re2_mean",
        "Re2_sd",
        "Re2_sem",
        "mean_free_path",
        "mean_free_path_sd",
        "mean_free_path_sem",
        "mean_free_path_fraction",
        "mean_free_path_fraction_sem",
        "seconds",
        "reliable",
    ]
    with open(args.out, "w", newline="") as handle:
        # Scrive un CSV piatto, facile da confrontare con i risultati pivot.
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.out}")
    if args.thesis:
        if (lengths != expected_lengths or args.walks != expected_walks
                or args.seed != expected_seed):
            raise SystemExit("--thesis arguments do not match the central protocol")
        refine_thesis_data(args.out, force=args.force_refinement,
                           out_dir=args.refinement_out,
                           budgets=refinement_budgets,
                           refinement_seed=refinement_seed)
    print(f"Total seconds: {time.perf_counter() - start_all:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
