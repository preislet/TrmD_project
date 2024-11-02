# Summary
This file will work as a log which will summarize team meetings

# Dates
- Deadline 1: 10.12.2024
- Midterm presentation 13.12.2024 10AM
- Deadline 2: 24.12.2024

## Meet 2 (1.11.2024)
- tasks were asigned to each team member:
    - Anna and Aneta - learn how to do target refinement via command line (possibly pyMol through python, removing waters, adding protonation.. preparing structures for further analysis), MSA of known target proteins and their variations 
    - Marek - classical virtual screening based on molecules from known databases
    - Tomáš - learn how to operate with ChEMBL and fetch approved molecules or, generally, our molecules of interest
    - Tom - virtual screening using generative approaches

## MeetEU (21.10.2024)
- Tomáš was assigned to read Evaluating druggability..
- Marián added Evolutionary repair, Development of inhibitors against Mycobacterium, Prospects of Indole, Unraveling allostery, Mg2+-dependent methyl transfer
- Tom's idea to dock right at the binding site of both AdoMet and tRNA with one ligand or double attack with two ligands
- Downside of bigger molecules is that they are expensive and they have many possible conformations (torsions...)
- Using covalent bond is riskant (toxicity, the ability to pass regulative pharmaceutical process)
- The problem of unideally resolved structural files (alignment of partially resolved and predicted structure) + ligand's end dynamics (AdoMet) and its possible bond/interaction within the structure (Marián thinks its transitional and does not stabilize the structure at all)
- Can the monomeric protein bind ligands?
- We want the structures of monomer, dimer, and all the other AlphaFold-predicted combinations (mono, di, holo, apo, and with ions - Mg) - basically get everything that's possible
    - files should be added in the shared discord server (we can then think of using MD to get other conformations/dynamic states)
    - review of all the described structures was assigned to Marián's student
- We found a file where the knot is structurally described (we have to look at its determined confidence; but unluckily it is a mutant in a position which can possibly interact with the knot itself)
    - we need to find out why was the mutation imposed
- We discussed a problem with big molecules in oficially-approved-drugs dataset

## Meet 1 (16.10.2024)
- [x] Random state: 42
- [ ] Keep documentation and comments on point
- [ ] Create workflow diagram
- [ ] How TrmD works + Mechanisms that we want to use
- [ ] Locate binding sites by p2rank
- [ ] Analysis of conformational changes in AHoJ-DB
- [ ] Protein-ligand corrections (SQM + ML → the best (PM6-ML))
- [ ] Approaches:
    - small molecules: destabilization of complex, work with Mg+2 (Ca+2, Mn+2), work with knot and AdoMet interaction,
    - proteins: destabilization of complex, mechanical block, PTM, work with C-terminal where tRNA binds (and is similar to Trp repressor),

- [ ] Settle on which approach to take or at least rank them
- [ ] Get relevant literature to present on Monday meeting (e.g. which supports our approach)
- [ ] All of us should take a look at the diagram (suggest corrections)
- [ ] As a first step gather databases to use in our workflow
