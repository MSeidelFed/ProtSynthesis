def count_gly_ser(peptide):

    peptide = peptide.split('_')

    if len(peptide) != 3:
        raise Exception('Parsing Error')
    else:
        peptide = peptide[1]

    modification = False
    gs_count = 0

    for aa in peptide:
        if aa == '[':
            modification = True
        elif aa == ']':
            modification = False
        elif modification == False:
            if aa == 'G' or aa =='S':
                gs_count += 1
        else: 
            continue

    return gs_count

        

