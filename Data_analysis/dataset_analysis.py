#!/usr/bin/env python
# coding: utf-8

## Packages
# import numpy as np
import pandas as pd
import math
import re
import argparse


# ================================== IMPORT DATASETS ==================================

def _datasets(file):  # srcfile/tgtfile
    data = pd.read_csv(file, header=None, names=['Reaction Smiles'], dtype=str)
    return data


def import_datasets(file_list):  # [train_srcf, train_tgtf, val_srcf, val_tgtf, test_srcf, test_tgtf]
    datasets_list = []
    for file in file_list:
        datasets_list.append(_datasets(file))
    print('''\n\nTRAIN\nSource and Target datasets loaded.......\n\n
             \nVALID\nSource and Target datasets loaded.......\n\n
             \nTEST\nSource and Target datasets loaded.......\n\n
             \n__________________________________________________________________''')
    return datasets_list  # [train_src, train_tgt, val_src, val_tgt, test_src, test_tgt]


# ================================== PROCESSING DATASETS ==================================

### Space removal
def _space_removal(dataset):
    dataframe = pd.DataFrame(dataset['Reaction Smiles'].str.replace(' ', ''))
    return dataframe


def space_removal(datasets_list):  # [train_src, train_tgt, val_src, val_tgt, test_src, test_tgt]
    for i in range(len(datasets_list)):
        datasets_list[i] = _space_removal(datasets_list[i])
    print('''\n\nTRAIN\nSource and Target datasets spaces removed.......\n\n
             \nVALID\nSource and Target datasets spaces removed.......\n\n
             \nTEST\nSource and Target datasets spaces removed.......\n\n
             \n__________________________________________________________________''')
    return datasets_list  # [train_src, train_tgt, val_src, val_tgt, test_src, test_tgt]


# ================================== SMILES ANALYSIS ==================================

## Auxiliary functions
def calc_state(total):
    total_1 = math.ceil(total / 100)
    total_125 = math.ceil(total / 8)
    total_25 = math.ceil(total / 4)
    total_375 = math.ceil(total * 3 / 8)
    total_50 = math.ceil(total / 2)
    total_625 = math.ceil(total * 5 / 8)
    total_75 = math.ceil(total * 3 / 4)
    total_875 = math.ceil(total * 7 / 8)
    return (total_1, total_125, total_25, total_375, total_50, total_625, total_75, total_875)


def state_report(total, total_1, total_125, total_25, total_375, total_50, total_625, total_75, total_875, index):
    '''Auxiliary function'''
    if index == total_1:
        print('= 1%')
    elif index == total_125:
        print('===== 12.5%')
    elif index == total_25:
        print('========== 25%')
    elif index == total_375:
        print('=============== 37.5%')
    elif index == total_50:
        print('==================== 50%')
    elif index == total_625:
        print('========================= 62.5%')
    elif index == total_75:
        print('============================== 75%')
    elif index == total_875:
        print('=================================== 87.5%')
    elif index == total:
        print('======================================== 100%\nDone')


def join_all(df_src, df_tgt):
    '''Join source and target smiles'''
    print("\n\nJoin Source and Target SMILES\nSource: ")
    print(df_src.shape)
    print("\nTarget: ")
    print(df_tgt.shape)
    df_reaction = df_src + '.' + df_tgt
    print("\nAll: ")
    print(df_reaction.shape)
    return df_reaction


def df_conversion(smiles_list):
    '''Conversion of smiles list into dataframe and drop duplicates'''
    df = pd.DataFrame(smiles_list, columns=['Smiles'])
    df = df.drop_duplicates()
    return df


# SMILES count
def smiles_count(df):
    '''Count number of smiles with and without stereochemical information'''
    counts = pd.DataFrame(columns=['N_smiles', 'N_stereo', 'N_nstereo'])
    smiles = []
    smiles_stereo = []
    smiles_nstereo = []
    # auxiliary variables
    total = len(df) - 1
    (total_1, total_125, total_25, total_375, total_50, total_625, total_75, total_875) = calc_state(total)
    print("\n\nInitializing...\n\nProcessing number of chiral molecules\n")
    # count number of smiles
    for index, row in df.iterrows():
        state_report(total, total_1, total_125, total_25, total_375, total_50, total_625, total_75, total_875, index)
        molecs = row['Reaction Smiles'].split('.')
        N_stereo = 0
        for mol in molecs:
            smiles.append(mol)
            if bool(re.search("\@", mol)):
                N_stereo += 1
                smiles_stereo.append(mol)
            else:
                smiles_nstereo.append(mol)
        counts.loc[index] = [len(molecs), N_stereo, len(molecs) - N_stereo]
    print('\nDropping duplicated SMILES')
    smiles = df_conversion(smiles)
    smiles_stereo = df_conversion(smiles_stereo)
    smiles_nstereo = df_conversion(smiles_nstereo)
    df = df.join(counts)
    return (df, smiles, smiles_stereo, smiles_nstereo)


#     return (counts, smiles, smiles_stereo, smiles_nstereo)

def smiles_all_concat(smiles_src, smiles_tgt):
    smiles_all = pd.concat([smiles_src, smiles_tgt], ignore_index=True)
    smiles_all = smiles_all.drop_duplicates()
    return smiles_all


def smiles_count_all(df_src, src_smiles, src_smiles_stereo, src_smiles_nstereo, df_tgt, tgt_smiles, tgt_smiles_stereo,
                     tgt_smiles_nstereo, df_reaction):
    print("\n\nInitializing...\n\nProcessing number of chiral molecules\n")
    state_report(1, 0, 0, 0, 0, 0, 0, 0, 0, 1)
    # adding counts
    cols = df_src.columns
    for col in cols[1:]:
        df_reaction[col] = df_src[col] + df_tgt[col]
        df_reaction['N_reagents'] = df_src['N_smiles']
        df_reaction['N_products'] = df_tgt['N_smiles']
    # SMILES
    print('\nDropping duplicated SMILES')
    all_smiles = smiles_all_concat(src_smiles, tgt_smiles)
    return (df_reaction, all_smiles)


# Reactions count
def reaction_loop(df, smile):
    n = 0
    smile = '.' + smile + '.'
    for reaction in df['Reaction Smiles']:
        reaction = '.' + reaction + '.'
        if smile in reaction:
            n += 1
    return n


def reaction_count(df_smiles, df_src, df_tgt, df_reaction):
    # count number of reactions per smile
    counts = pd.DataFrame(columns=['N_reactions', 'N_reagents', 'N_products'])
    # auxiliary variables
    total = len(df_smiles) - 1
    (total_1, total_125, total_25, total_375, total_50, total_625, total_75, total_875) = calc_state(total)
    print("\n\nInitializing...\n\nProcessing number of reactions per SMILES")
    # count number of reactions
    for index, smile in zip(range(len(df_smiles)), df_smiles['Smiles']):
        state_report(total, total_1, total_125, total_25, total_375, total_50, total_625, total_75, total_875, index)
        N_reactions = reaction_loop(df_reaction, smile)
        N_reagents = reaction_loop(df_src, smile)
        N_products = reaction_loop(df_tgt, smile)
        counts.loc[index] = [N_reactions, N_reagents, N_products]
    df_smiles = df_smiles.join(counts)
    return df_smiles


def verification(df_control, type_data, df_src, df_tgt, df_reaction, all_smiles, src_smiles_stereo, tgt_smiles_stereo,
                 src_smiles_nstereo, tgt_smiles_nstereo):
    '''Verification step regarding number of rows'''
    # DFs
    df_control.loc[len(df_control)] = [type_data, 'df_reaction', df_reaction.shape[0], df_reaction.shape[1],
                                       df_reaction.columns.values.tolist()]
    df_control.loc[len(df_control)] = [type_data, 'df_src', df_src.shape[0], df_src.shape[1],
                                       df_src.columns.values.tolist()]
    df_control.loc[len(df_control)] = [type_data, 'df_tgt', df_tgt.shape[0], df_tgt.shape[1],
                                       df_tgt.columns.values.tolist()]
    # SMILES
    # all
    df_control.loc[len(df_control)] = [type_data, 'all_smiles', all_smiles.shape[0], all_smiles.shape[1],
                                       all_smiles.columns.values.tolist()]
    # stereo
    df_control.loc[len(df_control)] = [type_data, 'src_smiles_stereo', src_smiles_stereo.shape[0],
                                       src_smiles_stereo.shape[1], src_smiles_stereo.columns.values.tolist()]
    df_control.loc[len(df_control)] = [type_data, 'tgt_smiles_stereo', tgt_smiles_stereo.shape[0],
                                       tgt_smiles_stereo.shape[1], tgt_smiles_stereo.columns.values.tolist()]
    # nstereo
    df_control.loc[len(df_control)] = [type_data, 'src_smiles_nstereo', src_smiles_nstereo.shape[0],
                                       src_smiles_nstereo.shape[1], src_smiles_nstereo.columns.values.tolist()]
    df_control.loc[len(df_control)] = [type_data, 'tgt_smiles_nstereo', tgt_smiles_nstereo.shape[0],
                                       tgt_smiles_nstereo.shape[1], tgt_smiles_nstereo.columns.values.tolist()]
    return df_control


# Counts
def df_count(df_src, df_tgt, df_control, data_type, file_name):
    print("\n\n__________________________" + data_type + "__________________________\n")
    df_reaction = join_all(df_src, df_tgt)
    # Source
    print("\n\n_ _ _ _ _ Source _ _ _ _ _")
    (df_src, src_smiles, src_smiles_stereo, src_smiles_nstereo) = smiles_count(df_src)
    # save
    df_src.to_csv(file_name + '_src.csv')
    src_smiles_stereo.to_csv(file_name + '_src_smiles_stereo.csv')
    src_smiles_nstereo.to_csv(file_name + '_src_smiles_nstereo.csv')
    # Target
    print("\n\n_ _ _ _ _ Target _ _ _ _ _")
    (df_tgt, tgt_smiles, tgt_smiles_stereo, tgt_smiles_nstereo) = smiles_count(df_tgt)
    # save
    df_tgt.to_csv(file_name + '_tgt.csv')
    tgt_smiles_stereo.to_csv(file_name + '_tgt_smiles_stereo.csv')
    tgt_smiles_nstereo.to_csv(file_name + '_tgt_smiles_nstereo.csv')
    # All
    print("\n\n_ _ _ _ _ All _ _ _ _ _")
    (df_reaction, all_smiles) = smiles_count_all(df_src, src_smiles, src_smiles_stereo, src_smiles_nstereo, df_tgt,
                                                 tgt_smiles, tgt_smiles_stereo, tgt_smiles_nstereo, df_reaction)
    # save
    df_reaction.to_csv(file_name + '_reaction.csv')
    # Smiles
    print("\n\n_ _ _ _ _ Smiles _ _ _ _ _")
    (all_smiles) = reaction_count(all_smiles, df_src, df_tgt, df_reaction)
    # save
    all_smiles.to_csv(file_name + '_all_smiles.csv')
    # verification step
    print("\n\n_ _ _ _ _ _ _ _ _ _ _\nControl step")
    df_control = verification(df_control, data_type, df_src, df_tgt, df_reaction, all_smiles, src_smiles_stereo,
                              tgt_smiles_stereo, src_smiles_nstereo, tgt_smiles_nstereo)
    # save
    print("\n\n_ _ _ _ _ _ _ _ _ _ _\n" + data_type + " dataframes successfully saved (folder: " + file_name[0:-6]+ "/)")
    print("\n________________________________________________________________________________________________________")
    return ([df_src, src_smiles_stereo, src_smiles_nstereo, df_tgt, tgt_smiles_stereo, tgt_smiles_nstereo, df_reaction,
             all_smiles], df_control)


def general_counts(datasets_list, dataset_type):  # train_src, train_tgt, val_src, val_tgt, test_src, test_tgt
    df_control = pd.DataFrame(columns=['Type data', 'DF name', 'N Rows', 'N Cols', 'Col name'])
    # type_of_data => src, src_smiles_stereo, src_smiles_nstereo, tgt, tgt_smiles_stereo, tgt_smiles_nstereo, reaction, all_smiles
    (train_list, df_control) = df_count(datasets_list[0], datasets_list[1], df_control, "TRAIN",
                                        dataset_type + '_analysis/train')
    (valid_list, df_control) = df_count(datasets_list[2], datasets_list[3], df_control, "VALID",
                                       dataset_type + '_analysis/test')
    (test_list, df_control) = df_count(datasets_list[4], datasets_list[5], df_control, "TEST",
                                        dataset_type + '_analysis/valid')
    df_control.to_csv(dataset_type + '_analysis/df_control.csv')
    print('\n\nSaving control information....\n\n\n\t\t\t\t\tALL DATASETS PROCESSED\n\n\n')
    return (train_list, valid_list, test_list, df_control)


def main(args):
    # file list - train_srcf, train_tgtf, val_srcf, val_tgtf, test_srcf, test_tgtf]
    # list datasets - train_src, train_tgt, val_src, val_tgt, test_src, test_tgt
    datasets_list = import_datasets(args.file_list)
    datasets_list = space_removal(datasets_list)
    (train_list, valid_list, test_list, df_control) = general_counts(datasets_list, args.data_type[0])


# fazer como argumento o tipo de dataset -STEREO etc - correr dentro da pasta dataset analysis
# dataset_type - STEREO_mixed
# file_list - [train_srcf, train_tgtf, val_srcf, val_tgtf, test_srcf, test_tgtf]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Process reaction smiles datasets, both src and tgt - train, test and validation datasets.')
    parser.add_argument('data_type', metavar='D', type=str, nargs=1, help='Type of dataset (e.g. STEREO_mixed_augm)')
    parser.add_argument('file_list', metavar='F', type=str, nargs='+', help='Path to dataset file/s.')
    args = parser.parse_args()
    main(args)
