# Leave credit at the top here because at some point in development I made a dumb mistake and deleted 
# a non-backed up version of this python file... Oof
# uncompyle6 version 3.8.0
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.18 |Anaconda, Inc.| (default, Apr 23 2020, 17:44:47) 
# [GCC 4.2.1 Compatible Clang 4.0.1 (tags/RELEASE_401/final)]

# Embedded file name: tcr_structParse.py
# Compiled at: 2022-05-02 12:35:38
import numpy as np, matplotlib.pyplot as pl, matplotlib as mpl, mdtraj as md, os
from Bio import pairwise2
import pandas, itertools

def convert_3Let(inp):
    first = True
    three_let = ['ALA', 'GLY', 'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'MET', 'CYS', 'PHE', 'TYR', 'THR', 'TRP', 'PRO', 'SER', 'LEU', 'VAL', 'HIS', 'ILE']
    sin_let = ['A', 'G', 'R', 'K', 'D', 'E', 'N', 'Q', 'M', 'C', 'F', 'Y', 'T', 'W', 'P', 'S', 'L', 'V', 'H', 'I']
    sin_final = []
    for i in inp:
        hold = []
        for scan in np.arange(len(three_let)):
            if i.lower() == three_let[scan].lower():
                hold = sin_let[scan]
                break

        if len(hold) == 0:
            continue
        if first:
            sin_final = hold
            first = False
        else:
            sin_final = np.hstack((sin_final, hold))

    if len(sin_final) == 0:
        return ()
    return sin_final

# Have HLA-DP as a standard
def get_chains(struct, mhc_class = 1,mhcID='HLA-DP'):
    trav_genes = pandas.read_csv('trav_human_full.csv')
    trbv_genes = pandas.read_csv('trbv_human_full.csv')
    if mhc_class == 1:
        mhc_seq = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWE'
    elif mhc_class == 2:
        if mhcID == 'HLA-DP':
            mhc_alpha = 'KADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHT'
        elif mhcID == 'HLA-DQ':
            mhc_alpha = 'ADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNST'
        elif mhcID == 'HLA-DR':
            mhc_alpha = 'KEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYT'
        mhc_beta = 'MVLLTSVVQGRATPENYVYQGRQECYAFNGTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEDYWNSQKDLLEEKRAVPDRVCRHNYELDEAVTLQ'

    table, bonds = struct.topology.to_dataframe()
    totChain = table['chainID'].drop_duplicates().values[(-1)]
    alpha_chain = -1
    beta_chain = -1
    if mhc_class == 1:
        # Due to a bunch of really jacked up PDBs, we need to build in a way to pull out distances when we aren't
        # selecting the proper chains of the asymmetric units (a bunch of "no contact" PDBs are due to this)
        mhc_chainList = []
        mhc_chain = -1
    elif mhc_class == 2:
        mhc_chainListA = []
        mhc_alpha_chain = -1
        mhc_chainListB = []
        mhc_beta_chain = -1
    alphaMax = 0
    betaMax = 0
    for chainNum in np.arange(totChain + 1):
        alpha_score = 0
        beta_score = 0
        chain_seq = [ residue for residue in struct.topology.chain(chainNum).residues ]
        seq_only = [ str(i)[:3] for i in chain_seq ]
        sinLet_seq = convert_3Let(seq_only)
        final_sequence = ('').join(sinLet_seq)
        if mhc_class == 1:
            aa = pairwise2.align.globalms(mhc_seq, final_sequence, 0.5, -0.1, -5, -0.5)
        elif mhc_class == 2:
            # Need to change the scoring a bit for class II, I think because of shorter sequence
            # in each individual matching
            aa_alpha = pairwise2.align.globalms(mhc_alpha, final_sequence, 1.5, -0.1, -5, -0.5)
            aa_beta = pairwise2.align.globalms(mhc_beta, final_sequence, 1.5, -0.1, -5, -0.5)
            if aa_alpha == [] or aa_beta == []:
                aa = []
            else:
                if aa_alpha[0][2] > aa_beta[0][2]:
                    aa = aa_alpha
                    classII_chain = 'alpha'
                else:
                    aa = aa_beta
                    classII_chain = 'beta'
        if aa == []:
            score = 0
        else:
            score = aa[0][2]
        if score > 0:
            if mhc_class == 1:
                mhc_chain = chainNum
                mhc_chainList = mhc_chainList + [chainNum]
            else:
                if classII_chain == 'alpha':
                    mhc_alpha_chain = chainNum
                    mhc_chainListA = mhc_chainListA + [chainNum]
                elif classII_chain == 'beta':
                    mhc_beta_chain = chainNum
                    mhc_chainListB = mhc_chainListB + [chainNum]
            continue
        for tcrA in trav_genes.values:
            temp_name = tcrA[0]
            temp_seq = tcrA[1]
            aa = pairwise2.align.localms(temp_seq, final_sequence, 0.5, -0.1, -5, -0.5)
            if aa == []:
                score = 0
            else:
                score = aa[0][2]
            if score > alpha_score:
                alpha_score = score
            if score > alphaMax:
                alpha_nameF = temp_name
                alphaMax = score

        for tcrB in trbv_genes.values:
            temp_name = tcrB[0]
            temp_seq = tcrB[1]
            aa = pairwise2.align.localms(temp_seq, final_sequence, 0.5, -0.1, -5, -0.5)
            if aa == []:
                score = 0
            else:
                score = aa[0][2]
            if score > beta_score:
                beta_score = score
            if score > betaMax:
                beta_nameF = temp_name
                betaMax = score
        # We're only going to return ONE TCR. This has worked well thus far
        if alpha_score > 20 and alpha_score >= alphaMax:
            alpha_chain = chainNum
        if beta_score > 20 and beta_score >= betaMax:
            beta_chain = chainNum
        # ALRIGHT SO A BIG CHANGE HERE. BECAUSE THE NEWER METHOD IS MUCH FASTER
        # LETS NOT BOTHER "BREAKING". LETS JUST MAKE A LIST OF EVERY MHC and EVERY TCR CHAIN

    if alpha_chain == -1:
        print('Cannot find TRAV match!')
        return ()
    if beta_chain == -1:
        print('Cannot find TRBV match!')
        return ()
    if mhc_class == 1:
        if mhc_chain == -1:
            print('Cannot find MHC match!')
            return ()
        return (alpha_chain, alpha_nameF, beta_chain, beta_nameF, mhc_chain,mhc_chainList)
    elif mhc_class == 2:
        if mhc_alpha_chain == -1:
            print('Cannot find MHCalpha match!')
            return ()
        if mhc_beta_chain == -1:
            print('Cannot find MHCbeta match!')
        return (alpha_chain, alpha_nameF, beta_chain, beta_nameF, mhc_alpha_chain,mhc_beta_chain,mhc_chainListA,mhc_chainListB)

def calc_process_classIdist(struct, tcr_chain, mhc_chain, alpha_nameF, beta_nameF,table,mhc_groupsel='sidechain',tcr_groupsel='sidechain',
                         ab='alpha', dist_cutoff=0.35,period=False):
    # New section of script to try to see if a new method would work
    # Define the MHC refDF
    mhc_top = struct.topology.select('chainid == ' + str(mhc_chain))
    chain_sub = table[(table['chainID'] == mhc_chain)].reset_index(drop=True)
    struct_serial_df = pandas.DataFrame(mhc_top,columns=['struct_ser'])
    mhc_refDF = pandas.concat([chain_sub,struct_serial_df],axis=1)
    # Define the TCR refDF
    tcr_top = struct.topology.select('chainid == ' + str(tcr_chain))
    tcr_tab = table[(table['chainID'] == tcr_chain)].reset_index(drop=True)
    tcr_serial_df = pandas.DataFrame(tcr_top,columns=['struct_ser'])
    tcr_refDF = pandas.concat([tcr_tab,tcr_serial_df],axis=1)

    # Try to edit this to make sure our numbers match from the get-go
    alpha1 = [55,56,59,60,63,66,67,70,74,77,80]
    pep_contact = [3,5,7,20,22,24,43,57,61,64,65,68,71,72,75,78,79,82,93,95,97,112,114,141,145,150,154,157,169]
    alpha2 = [143,144,147,148,149,152,153,156,160,161,164,165,167,168]

    hlaA_0101 = 'SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETL'

    trav_cdrs = pandas.read_csv('trav_human_cdrs.csv')
    trbv_cdrs = pandas.read_csv('trbv_human_cdrs.csv')
    trav_12seq = trav_cdrs[(trav_cdrs['gene'] == alpha_nameF)][['cdr1', 'cdr2']].values[0]
    trbv_12seq = trbv_cdrs[(trbv_cdrs['gene'] == beta_nameF)][['cdr1', 'cdr2']].values[0]
    tcr_sub = [ residue for residue in struct.topology.chain(tcr_chain).residues ]
    seq_only1 = [ str(i)[:3] for i in tcr_sub ]
    num_only1 = [ str(i)[3:] for i in tcr_sub ]
    sinLet_seq1 = convert_3Let(seq_only1)
    tcr_sequence = ('').join(sinLet_seq1)
    if ab == 'alpha':
        tcr_cdr1Start = tcr_sequence.find(trav_12seq[0])
        tcr_cdr2Start = tcr_sequence.find(trav_12seq[1])
        if tcr_cdr1Start == -1:
            tcr_cdr1Start = tcr_sequence.find(trav_12seq[0][0:3])
        if tcr_cdr1Start == -1:
            aa = pairwise2.align.localms(trav_12seq[0], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr1Start = aa[0][0].find(trav_12seq[0])
        if tcr_cdr2Start == -1:
            tcr_cdr2Start = tcr_sequence.find(trav_12seq[1][0:3])
        if tcr_cdr2Start == -1:
            aa = pairwise2.align.localms(trav_12seq[1], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr2Start = aa[0][0].find(trav_12seq[1])
        tcr_cdr1End = tcr_cdr1Start + len(trav_12seq[0])
        tcr_cdr2End = tcr_cdr2Start + len(trav_12seq[1])
    else:
        if ab == 'beta':
            tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0])
            tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1])
            if tcr_cdr1Start == -1:
                tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0][0:3])
            if tcr_cdr1Start == -1:
                aa = pairwise2.align.localms(trbv_12seq[0], tcr_sequence, 0.5, -0.1, -5, -0.5)
                tcr_cdr1Start = aa[0][0].find(trbv_12seq[0])
            if tcr_cdr2Start == -1:
                tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1][0:3])
            if tcr_cdr2Start == -1:
                aa = pairwise2.align.localms(trbv_12seq[1], tcr_sequence, 0.5, -0.1, -5, -0.5)
                tcr_cdr2Start = aa[0][0].find(trbv_12seq[1])
            tcr_cdr1End = tcr_cdr1Start + len(trbv_12seq[0])
            tcr_cdr2End = tcr_cdr2Start + len(trbv_12seq[1])

    mhc_sub = [ residue for residue in struct.topology.chain(mhc_chain).residues ]
    seq_only3 = [ str(i)[:3] for i in mhc_sub ]
    num_only3 = [ str(i)[3:] for i in mhc_sub ]
    sinLet_seq3 = convert_3Let(seq_only3)
    mhc_sequence = ('').join(sinLet_seq3)

    aa = pairwise2.align.localms(hlaA_0101,mhc_sequence,1.0, -0.1, -5,-0.5)
    num_shift = aa[0][0].find('SHSMRYF')
    resid_shift = 1 - int(num_only3[0])

    cdr1_len = tcr_cdr1End - tcr_cdr1Start
    cdr2_len = tcr_cdr2End - tcr_cdr2Start

    #############################################################################################################################
    # ALRIGHT SO WE REMOVED ALL OF THE REGISTER IDENTIFICATION STUFF. TRY TO REPLACE IT WITH A LOOKUP DATAFRAME
    # WHERE WE WERE CALLING FOR DIFFERENCES, INSTEAD CALL OUT THE EXACT MATCHES FURTHER DOWN
    #############################################################################################################################

    cdr1 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr1Start] + ' to ' + num_only1[tcr_cdr1End] + ') and '+tcr_groupsel)
    cdr2 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr2Start] + ' to ' + num_only1[tcr_cdr2End] + ') and '+tcr_groupsel)

    mhc_struct = struct.topology.select('chainid == ' + str(mhc_chain) + ' and '+mhc_groupsel)
    cdr1_mhc = list(itertools.product(cdr1, mhc_struct))
    cdr2_mhc = list(itertools.product(cdr2, mhc_struct))

    if len(cdr1) > 0:
        cdr1_dists = md.compute_distances(struct, atom_pairs=cdr1_mhc, periodic=period)
    else:
        print('BadSel CDR1'+ab)
    if len(cdr2) > 0:
        cdr2_dists = md.compute_distances(struct, atom_pairs=cdr2_mhc, periodic=period)
    else:
        print('BadSel CDR2'+ab)
    min_dist = dist_cutoff
    cdr1_pos = []
    cdr1_seldist = []

    if len(cdr1) > 0:
        for i in np.arange(len(cdr1_dists[0])):
            if cdr1_dists[0][i] < min_dist:
                cdr1_pos = cdr1_pos + [i]
                cdr1_seldist = cdr1_seldist + [cdr1_dists[0][i]]
    
    cdr2_pos = []
    cdr2_seldist = []
    if len(cdr2) > 0:
        for i in np.arange(len(cdr2_dists[0])):
            if cdr2_dists[0][i] < min_dist:
                cdr2_pos = cdr2_pos + [i]
                cdr2_seldist = cdr2_seldist + [cdr2_dists[0][i]]

    #################################################################################################
    #################################################################################################
    first = True
    for k in np.arange(len(cdr1_pos)):
        pairs = cdr1_pos[k]
        tcr_index = cdr1_mhc[pairs][0]
        mhc_index = cdr1_mhc[pairs][1]
        # Again hopefully this should work...
        tcr_ID = tcr_refDF[(tcr_refDF['struct_ser'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = mhc_refDF[(mhc_refDF['struct_ser'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]
        # Yea check for multiple hits just in case
        if np.shape(np.shape(tcr_ID))[0] != 1:
            print('Mutliple indices for one atom')
        if np.shape(np.shape(mhc_ID))[0] != 1:
            print('Mutliple indices for one atom')

        # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(num_only3)):
            if int(num_only3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhc_sequence[start_seq:start_seq+5]
        if hold_seq[0] != convert_3Let([mhc_ID[0]]):
            print('whoops')
        seq_loc = aa[0][1].find(hold_seq)

        for num in np.arange(len(alpha1)):
            if alpha1[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha1[num]

        for num in np.arange(len(alpha2)):
            if alpha2[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha2[num]
        ##################################

        if first:
            pre_df = np.hstack((tcr_ID, mhc_ID, cdr1_seldist[k], 'cdr1' + ab[0]))
            first = False
        else:
            pre_df = np.vstack((pre_df, np.hstack((tcr_ID, mhc_ID, cdr1_seldist[k], 'cdr1' + ab[0]))))
    #################################################################################################
    #################################################################################################

    for k in np.arange(len(cdr2_pos)):
        pairs = cdr2_pos[k]
        tcr_index = cdr2_mhc[pairs][0]
        mhc_index = cdr2_mhc[pairs][1]
        # This is really the only difference between old and new script.
        # Just find the "struct ser" and go with it... Might need a check for multiple hits, if there are any?
        tcr_ID = tcr_refDF[(tcr_refDF['struct_ser'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = mhc_refDF[(mhc_refDF['struct_ser'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]
        # Yea check for multiple hits just in case
        if np.shape(np.shape(tcr_ID))[0] != 1:
            print('Mutliple indices for one atom')
        if np.shape(np.shape(mhc_ID))[0] != 1:
            print('Mutliple indices for one atom')

        # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(num_only3)):
            if int(num_only3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhc_sequence[start_seq:start_seq+5]
        if hold_seq[0] != convert_3Let([mhc_ID[0]]):
            print('whoops')
        seq_loc = aa[0][1].find(hold_seq)

        for num in np.arange(len(alpha1)):
            if alpha1[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha1[num]

        for num in np.arange(len(alpha2)):
            if alpha2[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha2[num]
        ##################################

        if first:
            pre_df = np.hstack((tcr_ID, mhc_ID, cdr2_seldist[k], 'cdr2' + ab[0]))
            first = False
        else:
            pre_df = np.vstack((pre_df, np.hstack((tcr_ID, mhc_ID, cdr2_seldist[k], 'cdr2' + ab[0]))))

    if first == True:
        #print("No contacts found")
        return()
    final_df = pandas.DataFrame(pre_df)
    if np.shape(final_df)[1] != 8:
        final_df = np.transpose(final_df)
    final_df.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop']]
    return(final_df)
    # okay decompiling tcr_structParse.pyc

def calc_process_classIIdist(struct, tcr_chain, mhc_alpha_chain, mhc_beta_chain, alpha_nameF, beta_nameF,
table, ab='alpha', dist_cutoff=0.35,mhcID = 'HLA-DP',period=False,mhc_groupsel='sidechain',tcr_groupsel='sidechain'):
    # New section to try to get our new dataframe lookup (faster and more accurate)
    # distance processor that already works well for class I
    mhc_topA = struct.topology.select('chainid == ' + str(mhc_alpha_chain))
    structA_serial_df = pandas.DataFrame(mhc_topA,columns=['struct_serA'])
    mhc_topB = struct.topology.select('chainid == ' + str(mhc_beta_chain))
    structB_serial_df = pandas.DataFrame(mhc_topB,columns=['struct_serB'])
    chainA_sub = table[(table['chainID'] == mhc_alpha_chain)].reset_index(drop=True)
    chainB_sub = table[(table['chainID'] == mhc_beta_chain)].reset_index(drop=True)
    mhcA_refDF = pandas.concat([chainA_sub,structA_serial_df],axis=1)
    mhcB_refDF = pandas.concat([chainB_sub,structB_serial_df],axis=1)
    # Define the TCR refDF
    tcr_top = struct.topology.select('chainid == ' + str(tcr_chain))
    tcr_tab = table[(table['chainID'] == tcr_chain)].reset_index(drop=True)
    tcr_serial_df = pandas.DataFrame(tcr_top,columns=['struct_ser'])
    tcr_refDF = pandas.concat([tcr_tab,tcr_serial_df],axis=1)

    trav_cdrs = pandas.read_csv('trav_human_cdrs.csv')
    trbv_cdrs = pandas.read_csv('trbv_human_cdrs.csv')
    trav_12seq = trav_cdrs[(trav_cdrs['gene'] == alpha_nameF)][['cdr1', 'cdr2']].values[0]
    trbv_12seq = trbv_cdrs[(trbv_cdrs['gene'] == beta_nameF)][['cdr1', 'cdr2']].values[0]
    tcr_sub = [ residue for residue in struct.topology.chain(tcr_chain).residues ]
    seq_only1 = [ str(i)[:3] for i in tcr_sub ]
    num_only1 = [ str(i)[3:] for i in tcr_sub ]
    sinLet_seq1 = convert_3Let(seq_only1)
    tcr_sequence = ('').join(sinLet_seq1)
    if ab == 'alpha':
        tcr_cdr1Start = tcr_sequence.find(trav_12seq[0])
        tcr_cdr2Start = tcr_sequence.find(trav_12seq[1])
        if tcr_cdr1Start == -1:
            tcr_cdr1Start = tcr_sequence.find(trav_12seq[0][0:3])
        if tcr_cdr1Start == -1:
            aa = pairwise2.align.localms(trav_12seq[0], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr1Start = aa[0][0].find(trav_12seq[0])
        if tcr_cdr2Start == -1:
            tcr_cdr2Start = tcr_sequence.find(trav_12seq[1][0:3])
        if tcr_cdr2Start == -1:
            aa = pairwise2.align.localms(trav_12seq[1], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr2Start = aa[0][0].find(trav_12seq[1])
        tcr_cdr1End = tcr_cdr1Start + len(trav_12seq[0])
        tcr_cdr2End = tcr_cdr2Start + len(trav_12seq[1])
    else:
        if ab == 'beta':
            tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0])
            tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1])
            if tcr_cdr1Start == -1:
                tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0][0:3])
            if tcr_cdr1Start == -1:
                aa = pairwise2.align.localms(trbv_12seq[0], tcr_sequence, 0.5, -0.1, -5, -0.5)
                tcr_cdr1Start = aa[0][0].find(trbv_12seq[0])
            if tcr_cdr2Start == -1:
                tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1][0:3])
            if tcr_cdr2Start == -1:
                aa = pairwise2.align.localms(trbv_12seq[1], tcr_sequence, 0.5, -0.1, -5, -0.5)
                tcr_cdr2Start = aa[0][0].find(trbv_12seq[1])
            tcr_cdr1End = tcr_cdr1Start + len(trbv_12seq[0])
            tcr_cdr2End = tcr_cdr2Start + len(trbv_12seq[1])

    mhc_alpha_sub = [ residue for residue in struct.topology.chain(mhc_alpha_chain).residues ]
    mhc_beta_sub = [ residue for residue in struct.topology.chain(mhc_beta_chain).residues ]
    mhcalpha_seqonly3 = [ str(i)[:3] for i in mhc_alpha_sub ];mhcbeta_seqonly3 = [ str(i)[:3] for i in mhc_beta_sub ]
    mhcalpha_numonly3 = [ str(i)[3:] for i in mhc_alpha_sub ];mhcbeta_numonly3 = [ str(i)[3:] for i in mhc_beta_sub ]
    alphasinLet_seq3 = convert_3Let(mhcalpha_seqonly3); betasinLet_seq3 = convert_3Let(mhcbeta_seqonly3)
    mhcalpha_sequence = ('').join(alphasinLet_seq3); mhcbeta_sequence = ('').join(betasinLet_seq3)

    if mhcID == 'HLA-DP':
        mhc_alpha = 'KADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHT'
    elif mhcID == 'HLA-DQ':
        mhc_alpha = 'ADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNST'
    elif mhcID == 'HLA-DR':
        mhc_alpha = 'KEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYT'
    mhc_beta = 'MVLLTSVVQGRATPENYVYQGRQECYAFNGTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEDYWNSQKDLLEEKRAVPDRVCRHNYELDEAVTLQ'

    alpha = [51, 53, 55, 56, 58, 59, 62, 63, 66, 69, 70, 73, 74]
    beta = [63, 66, 67, 70, 71, 73, 74, 76, 77, 80, 83, 84, 87, 88, 91, 92]

    aa = pairwise2.align.localms(mhc_alpha,mhcalpha_sequence,1.0, -0.1, -5,-0.5)
    num_shift_alpha = aa[0][0].find(mhc_alpha[0:5])
    resid_shift_alpha = 1 - int(mhcalpha_numonly3[0])

    bb = pairwise2.align.localms(mhc_beta,mhcbeta_sequence,1.0, -0.1, -5,-0.5)
    num_shift_beta = bb[0][0].find(mhc_beta[0:5])
    resid_shift_beta = 1 - int(mhcbeta_numonly3[0])
    
    cdr1_len = tcr_cdr1End - tcr_cdr1Start
    cdr2_len = tcr_cdr2End - tcr_cdr2Start
    
    #############################################################################################################################
    # ALRIGHT SO WE REMOVED ALL OF THE REGISTER IDENTIFICATION STUFF. TRY TO REPLACE IT WITH A LOOKUP DATAFRAME
    # WHERE WE WERE CALLING FOR DIFFERENCES, INSTEAD CALL OUT THE EXACT MATCHES FURTHER DOWN
    #############################################################################################################################

    cdr1 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr1Start] + ' to ' + num_only1[tcr_cdr1End] + ') and '+tcr_groupsel)
    cdr2 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr2Start] + ' to ' + num_only1[tcr_cdr2End] + ') and '+tcr_groupsel)

    mhc_alpha_struct = struct.topology.select('chainid == ' + str(mhc_alpha_chain) + ' and '+mhc_groupsel)
    mhc_beta_struct = struct.topology.select('chainid == ' + str(mhc_beta_chain) + ' and '+mhc_groupsel)
    cdr1_mhcalpha = list(itertools.product(cdr1, mhc_alpha_struct))
    cdr1_mhcbeta = list(itertools.product(cdr1, mhc_beta_struct))
    cdr2_mhcalpha = list(itertools.product(cdr2, mhc_alpha_struct))
    cdr2_mhcbeta = list(itertools.product(cdr2, mhc_beta_struct))
    
    if len(cdr1) > 0:
        cdr1_dists_alpha = md.compute_distances(struct, atom_pairs=cdr1_mhcalpha, periodic=period)
        cdr1_dists_beta = md.compute_distances(struct, atom_pairs=cdr1_mhcbeta, periodic=period)
    else:
        return 'BadSel'
    if len(cdr2) > 0:
        cdr2_dists_alpha = md.compute_distances(struct, atom_pairs=cdr2_mhcalpha, periodic=period)
        cdr2_dists_beta = md.compute_distances(struct, atom_pairs=cdr2_mhcbeta, periodic=period)
    else:
        return 'BadSel'

    min_dist = dist_cutoff
    cdr1_pos_alpha = []; cdr1_seldist_alpha = []
    cdr1_pos_beta = []; cdr1_seldist_beta = []

    if len(cdr1) > 0:
        for i in np.arange(len(cdr1_dists_alpha[0])):
            if cdr1_dists_alpha[0][i] < min_dist:
                cdr1_pos_alpha = cdr1_pos_alpha + [i]
                cdr1_seldist_alpha = cdr1_seldist_alpha + [cdr1_dists_alpha[0][i]]

        for i in np.arange(len(cdr1_dists_beta[0])):
            if cdr1_dists_beta[0][i] < min_dist:
                cdr1_pos_beta = cdr1_pos_beta + [i]
                cdr1_seldist_beta = cdr1_seldist_beta + [cdr1_dists_beta[0][i]]

    cdr2_pos_alpha = []; cdr2_seldist_alpha = []
    cdr2_pos_beta = []; cdr2_seldist_beta = []

    if len(cdr2) > 0:
        for i in np.arange(len(cdr2_dists_alpha[0])):
            if cdr2_dists_alpha[0][i] < min_dist:
                cdr2_pos_alpha = cdr2_pos_alpha + [i]
                cdr2_seldist_alpha = cdr2_seldist_alpha + [cdr2_dists_alpha[0][i]]

        for i in np.arange(len(cdr2_dists_beta[0])):
            if cdr2_dists_beta[0][i] < min_dist:
                cdr2_pos_beta = cdr2_pos_beta + [i]
                cdr2_seldist_beta = cdr2_seldist_beta + [cdr2_dists_beta[0][i]]

    # The big changes are in this section, but really we're just changing how we find the TCRID
    # And the MHCID. Everything below that stays the same
    #################################################################################################
    #################################################################################################
    first_alpha = True
    for k in np.arange(len(cdr1_pos_alpha)):
        pairs = cdr1_pos_alpha[k]
        tcr_index = cdr1_mhcalpha[pairs][0]
        mhc_index = cdr1_mhcalpha[pairs][1]
        # These next 7 lines are the only new bits
        tcr_ID = tcr_refDF[(tcr_refDF['struct_ser'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = mhcA_refDF[(mhcA_refDF['struct_serA'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]
        # Yea check for multiple hits just in case
        if np.shape(np.shape(tcr_ID))[0] != 1:
            print('Mutliple indices for one atom')
        if np.shape(np.shape(mhc_ID))[0] != 1:
            print('Mutliple indices for one atom')

          # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcalpha_numonly3)):
            if int(mhcalpha_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcalpha_sequence[start_seq:start_seq+5]
        findit = aa[0][1].find(hold_seq)
        find_alpha=aa[0][0][findit:findit+5]
        for num in alpha:
            stand_alpha=mhc_alpha[num:num+5]
            if find_alpha == stand_alpha:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        ##################################
        if mhc_ID[1] == 'DROP':
            continue

        if first_alpha:
            pre_df_alpha = np.hstack((tcr_ID, mhc_ID, cdr1_seldist_alpha[k], 'cdr1' + ab[0]))
            first_alpha = False
        else:
            pre_df_alpha = np.vstack((pre_df_alpha, np.hstack((tcr_ID, mhc_ID, cdr1_seldist_alpha[k], 'cdr1' + ab[0]))))
    
    first_beta = True
    for k in np.arange(len(cdr1_pos_beta)):
        pairs = cdr1_pos_beta[k]
        tcr_index = cdr1_mhcbeta[pairs][0]
        mhc_index = cdr1_mhcbeta[pairs][1]
        tcr_ID = tcr_refDF[(tcr_refDF['struct_ser'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = mhcB_refDF[(mhcB_refDF['struct_serB'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]
        # Yea check for multiple hits just in case
        if np.shape(np.shape(tcr_ID))[0] != 1:
            print('Mutliple indices for one atom')
        if np.shape(np.shape(mhc_ID))[0] != 1:
            print('Mutliple indices for one atom')

         # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcbeta_numonly3)):
            if int(mhcbeta_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcbeta_sequence[start_seq:start_seq+5]
        findit = bb[0][1].find(hold_seq)
        find_beta=bb[0][0][findit:findit+5]
        for num in beta:
            stand_beta=mhc_beta[num:num+5]
            if find_beta == stand_beta:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        if mhc_ID[1] == 'DROP':
            continue
        ##################################

        if first_beta:
            pre_df_beta = np.hstack((tcr_ID, mhc_ID, cdr1_seldist_beta[k], 'cdr1' + ab[0]))
            first_beta = False
        else:
            pre_df_beta = np.vstack((pre_df_beta, np.hstack((tcr_ID, mhc_ID, cdr1_seldist_beta[k], 'cdr1' + ab[0]))))

    for k in np.arange(len(cdr2_pos_alpha)):
        pairs = cdr2_pos_alpha[k]
        tcr_index = cdr2_mhcalpha[pairs][0]
        mhc_index = cdr2_mhcalpha[pairs][1]
        tcr_ID = tcr_refDF[(tcr_refDF['struct_ser'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = mhcA_refDF[(mhcA_refDF['struct_serA'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]
        # Yea check for multiple hits just in case
        if np.shape(np.shape(tcr_ID))[0] != 1:
            print('Mutliple indices for one atom')
        if np.shape(np.shape(mhc_ID))[0] != 1:
            print('Mutliple indices for one atom')

         # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcalpha_numonly3)):
            if int(mhcalpha_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcalpha_sequence[start_seq:start_seq+5]
        findit = aa[0][1].find(hold_seq)
        find_alpha=aa[0][0][findit:findit+5]
        for num in alpha:
            stand_alpha=mhc_alpha[num:num+5]
            if find_alpha == stand_alpha:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        ##################################
        if mhc_ID[1] == 'DROP':
            continue

        if first_alpha:
            pre_df_alpha = np.hstack((tcr_ID, mhc_ID, cdr2_seldist_alpha[k], 'cdr2' + ab[0]))
            first_alpha = False
        else:
            pre_df_alpha = np.vstack((pre_df_alpha, np.hstack((tcr_ID, mhc_ID, cdr2_seldist_alpha[k], 'cdr2' + ab[0]))))

    for k in np.arange(len(cdr2_pos_beta)):
        pairs = cdr2_pos_beta[k]
        tcr_index = cdr2_mhcbeta[pairs][0]
        mhc_index = cdr2_mhcbeta[pairs][1]
        tcr_ID = tcr_refDF[(tcr_refDF['struct_ser'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = mhcB_refDF[(mhcB_refDF['struct_serB'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]
        # Yea check for multiple hits just in case
        if np.shape(np.shape(tcr_ID))[0] != 1:
            print('Mutliple indices for one atom')
        if np.shape(np.shape(mhc_ID))[0] != 1:
            print('Mutliple indices for one atom')

         # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcbeta_numonly3)):
            if int(mhcbeta_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcbeta_sequence[start_seq:start_seq+5]
        findit = bb[0][1].find(hold_seq)
        find_beta=bb[0][0][findit:findit+5]
        for num in beta:
            stand_beta=mhc_beta[num:num+5]
            if find_beta == stand_beta:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        if mhc_ID[1] == 'DROP':
            continue
        ##################################

        if first_beta:
            pre_df_beta = np.hstack((tcr_ID, mhc_ID, cdr2_seldist_beta[k], 'cdr2' + ab[0]))
            first_beta = False
        else:
            pre_df_beta = np.vstack((pre_df_beta, np.hstack((tcr_ID, mhc_ID, cdr2_seldist_beta[k], 'cdr2' + ab[0]))))

    if first_alpha == True and first_beta == True:
        #print('no found contacts')
        return ()

    if first_alpha:
        final_df_alpha = []
    else:
        final_df_alpha = pandas.DataFrame(pre_df_alpha)
        if np.shape(final_df_alpha)[1] != 8:
            final_df_alpha = np.transpose(final_df_alpha)
        final_df_alpha.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop']]
    
    if first_beta:
        final_df_beta = []
    else:
        final_df_beta = pandas.DataFrame(pre_df_beta)
        if np.shape(final_df_beta)[1] != 8:
            final_df_beta = np.transpose(final_df_beta)
        final_df_beta.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop']]
    
    return(final_df_alpha,final_df_beta)

# This whole bit of code is to break down on a by-residue basis the AIMS Scoring for TCR-MHC Interactions
def get_germline_breakdown(SCORES,tcr_names,hla_names,mhc_type = 'classI',ScoreAlpha1=True,ScoreAlpha2 = False,
                           len_weight=False,score_weight=False,clash=False,BADclash = False):
    # In this, we want to take the scores and identify WHERE
    # The positive interactions are coming from.
    # This is for if you want a subset
    if mhc_type == 'classI':
        alpha1 = [55,56,59,60,63,66,67,70,74,77,80]
        alpha2 = [143,144,147,148,149,152,153,156,160,161,164,165,167,168]
    elif mhc_type == 'classII_alpha':
        II_alpha_contacts = [51, 53, 55, 56, 58, 59, 62, 63, 66, 69, 70, 73, 74]
        alpha1 = II_alpha_contacts
        alpha2 = []
    elif mhc_type == 'classII_beta':
        II_beta_contacts = [63, 66, 67, 70, 71, 73, 74, 76, 77, 80, 83, 84, 87, 88, 91, 92]
        alpha1 = II_beta_contacts
        alpha2 = []
    
    output_hist = np.zeros((len(tcr_names),len(hla_names)))
    output_df = pandas.DataFrame(output_hist)
    output_df.columns = [i[0] for i in hla_names]
    output_df.index = [i[0] for i in tcr_names]
    First = True
    for i in hla_names:
        for j in tcr_names:
            if First:
                fin_name = [i[0],j[0]]
                First = False
            else:
                fin_name = np.vstack((fin_name,[i[0],j[0]]))
    
    noMatch = []
    first = True
    for seq in np.arange(len(SCORES)):
        lochla = fin_name[seq][0]
        loctcr = fin_name[seq][1]
        test_score = SCORES[seq]
        len1, len2 = np.shape(test_score)
        save_coords = []
        for i in np.arange(len1):
            # Break when you are reaching past the end
            if i + 2 >= len1:
                break
            if ScoreAlpha1:
                #do it for the alpha1 loop
                for j in np.arange(len(alpha1)):
                    # Break when you are past the end
                    if j + 2 >= len(alpha1):
                        break
                    # Looking for contiguous regions
                    # of some positive interactions
                    # can change criteria later
                    test1 = test_score[i,alpha1[j]]
                    test2 = test_score[i+1,alpha1[j+1]]
                    test3 = test_score[i+2,alpha1[j+2]]
                    if clash:
                        tf1 = (test1 < 0)
                        tf2 = (test2 < 0)
                        tf3 = (test3 < 0)
                    elif BADclash:
                        # only test tf1
                        tf1 = (test1 == -2)
                        tf2 = True
                        tf3 = True
                    else:
                        tf1 = (test1 > 0)
                        tf2 = (test2 > 0)
                        tf3 = (test3 > 0)
                    if tf1 & tf2 & tf3:
                        # Save only the middle coords
                        save_coords = save_coords + [[i+1,alpha1[j],alpha1[j+1],alpha1[j+2]]]

                        # Add to a counter using this complicated line
                        if score_weight:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + (test1+test2+test3)/3
                        elif BADclash:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + test1
                        else:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + 1
            if ScoreAlpha2:
                # Do it again for the alpha2 loop
                for j in np.arange(len(alpha2)):
                    # Break when you are past the end
                    if j + 2 >= len(alpha2):
                        break
                    # Looking for contiguous regions
                    # of some positive interactions
                    # can change criteria later
                    test1 = test_score[i,alpha2[j]]
                    test2 = test_score[i+1,alpha2[j+1]]
                    test3 = test_score[i+2,alpha2[j+2]]
                    if clash:
                        tf1 = (test1 < 0)
                        tf2 = (test2 < 0)
                        tf3 = (test3 < 0)
                    elif BADclash:
                        # only test tf1
                        tf1 = (test1 == -2)
                        tf2 = True
                        tf3 = True
                    else:
                        tf1 = (test1 > 0)
                        tf2 = (test2 > 0)
                        tf3 = (test3 > 0)
                    if tf1 & tf2 & tf3:
                        # Save only the middle coords
                        save_coords = save_coords + [[i+1,alpha2[j],alpha2[j+1],alpha2[j+2]]]

                        # Add to the counter above
                        if score_weight:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + (test1+test2+test3)/3
                        elif BADclash:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + test1
                        else:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + 1
                            
        # Before we move on, if you want to, scale all these by the CDR lengths
        # Obviously based on this simple counting metric, longer CDRs will have
        # higher interaction scores by default
        if len_weight:
            # len1 should give the length of the CDR of interest
            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla]/len1
        if save_coords == []:
            # There were no contiguous matches for a bunch of these sequences!
            # Save the list of these sequences to be sure
            noMatch = noMatch + [seq]
            continue
            
        if first:
            final_coords = save_coords
            first = False
        else:
            final_coords = np.vstack((final_coords,save_coords))
    if first:
        return(1,2,3)
    else:
        return(output_df,final_coords,noMatch)