import dictionaries

def flexcalc(protein, r):
    '''Function that calculates the residues b-factors according to their neighbors'''             
    b_query =  None
    query = protein[r]
    if r == 0:
        b_query = dictionaries.b_factor[query]
    elif r == len(protein) - 1:
        b_query = dictionaries.b_factor[query]
    else:
        if protein[r - 1] in dictionaries.rigid:
            if protein[r + 1] in dictionaries.rigid:
                b_query = dictionaries.b_factor_rr[query]
            elif protein[r + 1] in dictionaries.flex:
                b_query = dictionaries.b_factor_rf[query]
        elif protein[r - 1] in dictionaries.flex:
            if protein[r + 1] in dictionaries.rigid:
                b_query = dictionaries.b_factor_rf[query]
            elif protein[r + 1] in dictionaries.flex:
                b_query = dictionaries.b_factor_ff[query]
    return b_query

protein = 'XGATQSFQSVGDLTPAEKDLIRSTWDQLMTHRTGFVADVFIRIFHNDPTAQRKFPQMAGLSPAELRTSRQMHAHAIRVSALMTTYIDEMDTEVLPELLATLTRTHDKNHVGKKNYDLFGKVLMEAIKAELGVGFTKQVHDAWAKTFAIVQGVLITKHAS'
print(flexcalc(protein, 1))