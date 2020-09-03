def compare_two_sets(up_1,up_2,dn_1,dn_2):

    up_up = len(up_1.intersection(up_2))
    dn_dn = len(dn_1.intersection(dn_2))
    dn_up = len(dn_1.intersection(up_2))
    up_dn = len(up_1.intersection(dn_2))

    union = len(up_1.union(up_2).union(dn_1).union(dn_2))
    return (up_up + dn_dn-dn_up-up_dn)/(union)
