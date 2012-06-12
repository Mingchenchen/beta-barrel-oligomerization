def revise(pdbid):
    # In order for selections to group properly,
    # the name of the group must have the right case in a group.sele
    # expression
    pdbid = pdbid.upper()
    
    # Remove the residues in "sele" from "pdbid.removed"
    cmd.select(pdbid + '.removed', pdbid + '.removed & ! sele')
    
    # Recoler the protein so I can verify that everything worked as expected
    cmd.color('green', pdbid + '.molecule')
    cmd.color('red', pdbid + '.removed')