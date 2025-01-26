# Summary
This file will work as a log which will summarize team meetings

# Dates
- Deadline 1: 10.12.2024
- Midterm presentation 13.12.2024 10AM
- Deadline 2: 24.12.2024


## MeetEU (8.12.2024)
- prepare our presentation (which molecules, where to dock, how and why)

## MeetEU (2.12.2024)
- suggested how our workflow could look
- revalidation after screening (of the top rated molecules; SQM; ADR prediction?)
- fill the form for 4EU+ trip to Milan

## MeetEU (25.11.2024)
- proS seems to replace the function of TrmD
- the is not a homologous protein to TrmD which would replace the function
- we should ideally target TrmD and proS at once (take a look at the alignment if there are possible regions to target)
- paper showcase where they specifically mutated at the binding site but mutations completely outside of the pocket had effect as well (maybe it restricts dimerization or it causes the protein to fold badly)
- two (known) possible functional analogues
- idea to block dimerization and the binding of AdoMet at once

## MeetEU (18.11.2024)
- Tom has had a presentation on molecular docking algorithms
- he suggested that a deep learning approach should be used and then after the ligand is fitted a physics-based approach should be used as well to introduce realistic conformation of the ligand in the pocket
- for number of reasons UniDock seems to be the one to use out of physics-based methods and SurfDock and uniMolDock of deep learning ones

## MeetEU (11.11.2024)
- Marian told us about a study which highlights the biological importance of methylations - after the knockout of methyl transferase gene, another protein replaced its function
- Tom was assigned to have a presentation on different molecular docking techniques (algorithms)
- we should find homologous sequences (with differet sequence identity in the binding sites) - however, these proteins should still be related to methylation
- the teams should agree to take different approaches (AdoMet or tRNA binding site)

## MeetEU (6.11.2024)
- knocking out or mutating the gene of TrmD or inhibiting the targets was met with no effect (study says)
- two known small molecules suggested as good inhibitors (or they have great binding affinity at least)

## MeetEU (4.11.2024)
- crystal structures were presented
- we should compare images of ligand inner-interactions of resolved crystal structures from PDBsum
- Tom has shown that the AdoMet pocket is very rigid (in dimers as well as in monomers - we are not sure if monomers weren't just extracted from dimeric crystal structure)
- we agreed on docking into structure 4YVG which has the best resolution and is monomeric (if we would have wanted to dock into dimeric structure then AlphaFold would render useful and relevant)
- it is not completely out of question to dock into second pocket (tRNA binding site of dimeric structure)

## Team-meet 2 (1.11.2024)
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

## Team-meet 1 (16.10.2024)
- [x] Random state: 42
- [ ] Keep documentation and comments on point
- [x] Create workflow diagram
- [ ] How TrmD works + Mechanisms that we want to use
- [x] Locate binding sites by p2rank
- [ ] Analysis of conformational changes in AHoJ-DB
- [ ] Protein-ligand corrections (SQM + ML → the best (PM6-ML))
- [ ] Approaches:
    - small molecules: destabilization of complex, work with Mg+2 (Ca+2, Mn+2), work with knot and AdoMet interaction,
    - proteins: destabilization of complex, mechanical block, PTM, work with C-terminal where tRNA binds (and is similar to Trp repressor),

- [ ] Settle on which approach to take or at least rank them
- [x] Get relevant literature to present on Monday meeting (e.g. which supports our approach)
- [x] All of us should take a look at the diagram (suggest corrections)
- [x] As a first step gather databases to use in our workflow
