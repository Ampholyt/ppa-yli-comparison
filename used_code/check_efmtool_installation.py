import efmtool
import numpy as np

if __name__ == "__main__":
    # Network structure:
    # => A, A <=> B, B <=> C, A => C, C <=>
    S = np.array(
        [
            [1, -1, 0, -1, 0],
            [0, 1, -1, 0, 0],
            [0, 0, 1, 1, -1],
        ]
    )
    rev = [0, 1, 1, 0, 1]
    reactions = ["r1", "r2", "r3", "r4", "r5"]
    metabolites = ["A", "B", "C"]
    efms = efmtool.calculate_efms(S, rev, reactions, metabolites)

    print(f"Found {efms.shape[1]} EFMs.")
