import gemmi
path = "bad_phase_test/8ejo/8ejo.pdb"
s = gemmi.read_structure(path)

ns = gemmi.Structure()
ns.cell = s.cell
nm = gemmi.Model(s[0].name)
for c in s[0]:
    nc = gemmi.Chain(c.name)
    for r in c: 
        i = gemmi.find_tabulated_residue(r.name)
        if i.is_nucleic_acid() or i.is_water(): 
            continue

        nc.add_residue(r)
    nm.add_chain(nc)
ns.add_model(nm)
ns.write_pdb("bad_phase_test/8ejo/8ejo_removed_na_water.pdb")

    

