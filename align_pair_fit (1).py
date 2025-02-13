from pymol import cmd

cmd.load("N1.xyz", "N1")
cmd.load("N2.xyz", "N2")
cmd.load("P1.xyz", "P1")
cmd.load("P2.xyz", "P2")

cmd.show("sticks")
cmd.hide("spheres")

#label all, name+resi

cmd.pair_fit(
    "P2 and 28/O",
    "N2 and 28/O"
)
cmd.pair_fit(
    "P1 and 28/O",
    "N2 and 28/O"
)
cmd.pair_fit(
    "N1 and 28/O",
    "N2 and 28/O"
)
cmd.pair_fit(
    "P2 and (1/C, 2/C, 7/H, 3/C, 28/O, 29/H, 4/C, 8/H, 5/C, 9/H, 6/C, 10/N, 11/C 12/H 13/H 14/H 15/C 16/H 17/H 18/H )",
    "N2 and (1/C, 2/C, 7/H, 3/C, 28/O, 29/H, 4/C, 8/H, 5/C, 9/H, 6/C, 10/N, 11/C 12/H 13/H 14/H 15/C 16/H 17/H 18/H )"
)
cmd.pair_fit(
    "P1 and (1/C, 2/C, 7/H, 3/C, 28/O, 29/H, 4/C, 8/H, 5/C, 9/H, 6/C, 10/N, 11/C 12/H 13/H 14/H 15/C 16/H 17/H 18/H )",
    "N2 and (1/C, 2/C, 7/H, 3/C, 28/O, 29/H, 4/C, 8/H, 5/C, 9/H, 6/C, 10/N, 11/C 12/H 13/H 14/H 15/C 16/H 17/H 18/H )"
)
cmd.pair_fit(
    "N1 and (1/C, 2/C, 7/H, 3/C, 28/O, 29/H, 4/C, 8/H, 5/C, 9/H, 6/C, 10/N, 11/C 12/H 13/H 14/H 15/C 16/H 17/H 18/H )",
    "N2 and (1/C, 2/C, 7/H, 3/C, 28/O, 29/H, 4/C, 8/H, 5/C, 9/H, 6/C, 10/N, 11/C 12/H 13/H 14/H 15/C 16/H 17/H 18/H )"
)

N1 aligned with P1
