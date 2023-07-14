import argparse
import pickle

import pandas as pd

''' 
    TODO: 2023-07-14 - Add automatic plots
''' 

####################################################################
## Global variables ################################################
####################################################################

col_class_eve = 'class'
col_class_dark = 'Similarity_label'
col_id_dark = 'Contigs_ID'

x_range_eve = list(range(4, 52))
x_range_dark = list(range(2, 51))

feat_eves = [
    'X15', 'X16', 'X17', 'X18', 'X19', 'X20', 'X21', 'X22', 'X23', 'X24',
    'X25', 'X26', 'X27', 'X28', 'X29', 'X30', 'X31', 'X32', 'X33', 'X34',
    'X35', 'X.15', 'X.16', 'X.17', 'X.18', 'X.19', 'X.20', 'X.21', 'X.22',
    'X.23', 'X.24', 'X.25', 'X.26', 'X.27', 'X.28', 'X.29', 'X.30', 'X.31',
    'X.32', 'X.33', 'X.34', 'X.35', 'dens15to18', 'dens20to22',
    'dens25to29', 'ratiosi_pi', 'ratio_si', 'dens18to35'
]

df_pipe: pd.DataFrame

# Command line arguments
path_input: str
path_output: str
path_classifier: str

is_verbose = False

####################################################################
## Main ############################################################
####################################################################

def parse_input():
    ''' 
        TODO: 2023-07-14 - ADD Description
    '''

    # Create args parser
    parser = argparse.ArgumentParser(
        description='Runs Viral × EVEs classification of small RNA profile represented contigs and export result as CSV',
        add_help=True,
        allow_abbrev=True,
        exit_on_error=True,
    )

    parser.add_argument('--input', help='Path to input')
    parser.add_argument('--output', help='Path to export resulting .csv output')
    parser.add_argument('--classifier', help='Path to pickle classifier model')
    parser.add_argument('--verbose', action='store_true', help='Enable is_verbose mode')

    # Set global values
    args = parser.parse_args()

    path_input = args.input
    path_output = args.output
    path_classifier = args.classifier
    is_verbose = args.verbose
    
    df_pipe = pd.read_table(path_input)

def get_classification_df() -> tuple:
    ''' 
        Transform data to classify

        TODO: 2023-06-17 - ADD Description
    '''

    ids = df_pipe[col_id_dark] # Contig ID's
    X_dark = df_pipe.iloc[:, x_range_dark] # Feature Values
    y_dark = df_pipe[col_class_dark] # Labels

    if (is_verbose):
        n_samples, n_features = X_dark.shape
        
        print('-----------------------------------')
        print(f"n_samples: {n_samples} × n_features: {n_features}")
        print(f"ID's -> ids.shape {ids.shape} ({ids.unique().shape[0]} are unique)")
        print(f'y.shape ({col_class_dark}): {y_dark.shape}')
        print(f'X.shape: {X_dark.shape}')

        print('-----------------------------------')
        feat_dark = X_dark.columns
        print(f'feat_dark: {feat_dark}')
        print(f'set(X_dark.dtypes): {set(X_dark.dtypes)}')

        print('-----------------------------------')
        classes = list(y_dark.unique())
        class_counts = list(y_dark.value_counts())
        # print(f'classes: {classes}')
        # print(f'class_counts: {class_counts}')
        print(y_dark.value_counts())

    feat_dark = X_dark.columns
    feat_common = [f for f in feat_eves if f in feat_dark]

    # Create new data frame
    df_classif = pd.DataFrame()
    df_classif[col_id_dark] = df_pipe[col_id_dark]
    df_classif[col_class_dark] = df_pipe[col_class_dark]

    # Add features with proper names
    feat_classif = []

    for i in range(15, 35 + 1):
        feat_sense = f'X{i}'
        feat_classif.append(feat_sense)
        df_classif[feat_sense] = X_dark[f'{i}']

    for i in range(15, 35 + 1):
        feat_anti_sense = f'X.{i}'
        feat_classif.append(feat_anti_sense)
        df_classif[feat_anti_sense] = X_dark[f'-{i}']

    df_classif[feat_common] = X_dark[feat_common].copy()
    feat_classif += feat_common.copy()

    if (is_verbose):
        print('-----------------------------------')
        # df_classif = df_classif.copy()

        n_match_feats = sum(feat_classif == feat_eves)
        print(f'feat_eves: {len(feat_eves)}') 
        print(f'feat_classif: {len(feat_classif)}')
        print(f'n_match_feats: {n_match_feats}')

    return df_classif, feat_classif

def __main__() -> None:
    ''' 
        TODO: 2023-06-17 - ADD Description
    '''

    parse_input()

    # Transform data to classify
    df_classif, feat_classif = get_classification_df()

    # Run classification
    with open(path_classifier, 'rb') as file:
        classifier = pickle.load(file)

    y_hat = classifier.predict(df_classif[feat_classif])
    df_classif[col_class_eve] = y_hat.copy()
    
    # Export csv
    df_classif.to_csv(path_output)


if __name__ == "__main__" : __main__()