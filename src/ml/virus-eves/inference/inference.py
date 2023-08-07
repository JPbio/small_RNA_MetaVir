import pickle

import pandas as pd

# n_samples: 80 × n_features: 49
# ID's -> ids.shape (80,) (80 are unique)
# y.shape (Similarity_label): (80,)
# X.shape: (80, 49)
# classes: ['viral', 'nonviral', 'nohit']
# class_counts: [42, 24, 14]
# viral       42
# nonviral    24
# nohit       14
# Name: Similarity_label, dtype: int64
# feat_eves: 48
# feat_dark2: 48
# n_match_feats: 48

# config_example = {
    
#     # Input
#     'input_size': 0,
#     'n_reads': 0,
    
#     # Processing
#     'n_contigs': 0,
#     'n_contigs_gt_200': 0,

#     # Classification
#     'classification': {
#         'pipe': {
#             'viral': 0,
#             'nonviral': 0,
#             'nohit': 0,
#         },
#         'classifier': {
#             'viral': 0,
#             'eve': 0,
#         },
#         'comparison': {
#             'viral_viral': 0,
#             'viral_eve': 0,
            
#             'nonviral_viral': 0,
#             'nonviral_eve': 0,
            
#             'nohit_viral': 0,
#             'nohit_eve': 0,
#         },
#     },

#     # Output
#     'output_size': 0, # Tamanho da pasta de saída [GB]

#     # Times
#     'time': {
#         'total': 0, # REVIEW: 2023-06-17 - Shouldn't this be calculated?
#     },
# }

def get_input_x_y_ids(
    df_dark: pd.DataFrame, x_range_dark: list,
    col_id_dark: str, col_class_dark: str
) -> tuple:
    ''' 
        Separate: X _(features)_ × Y _(classes)_ ×  ID's _(contig ID's)_
        
        TODO: 2023-06-17 - ADD Description
    ''' 
    
    ids = df_dark[col_id_dark] # Contig ID's
    X_dark = df_dark.iloc[:, x_range_dark] # Feature Values
    y_dark = df_dark[col_class_dark] # Labels
    return X_dark, y_dark, ids

def get_classification_df(
    X_dark: pd.DataFrame, df_dark: pd.DataFrame, feat_eves: list,
    col_class_dark: str, col_id_dark: str
) -> tuple:
    ''' 
        Transform data to classify

        TODO: 2023-06-17 - ADD Description
    '''

    feat_dark = X_dark.columns
    feat_common = [f for f in feat_eves if f in feat_dark]

    # Create new data frame
    df_dark2 = pd.DataFrame()

    df_dark2[col_id_dark] = df_dark[col_id_dark]
    df_dark2[col_class_dark] = df_dark[col_class_dark]

    # Add features with proper names
    feat_dark2 = []

    for i in range(15, 35 + 1):
        feat_sense = f'X{i}'
        feat_dark2.append(feat_sense)
        df_dark2[feat_sense] = X_dark[f'{i}']

    for i in range(15, 35 + 1):
        feat_anti_sense = f'X.{i}'
        feat_dark2.append(feat_anti_sense)
        df_dark2[feat_anti_sense] = X_dark[f'-{i}']

    df_dark2[feat_common] = X_dark[feat_common].copy()
    feat_dark2 += feat_common.copy()

    return df_dark2, feat_dark2

# def classify(df_dark2: pd.DataFrame, feat_dark2: list, classifier) -> pd.DataFrame:
def classify(df_dark2: pd.DataFrame, feat_dark2: list, path_classifier: str) -> pd.DataFrame:
    ''' 
        Classify Dark Matter DB sequences
    ''' 

    with open(path_classifier, 'rb') as file:
        classifier = pickle.load(file)

    y_hat = classifier.predict(df_dark2[feat_dark2])
    return y_hat

def build_report(
    # df_dark2: pd.DataFrame, y_hat: pd.Series,
    df_dark2: pd.DataFrame,
    col_class_dark: str, col_class_eve: str, col_id_dark: str,
) -> pd.DataFrame:
    ''' 
        TODO: 2023-06-17 - ADD Description
    ''' 

    # df_dark2[col_class_eve] = y_hat.copy()
    # report_name = f'{lib_id}-report'

    cols_report = [col_id_dark, col_class_dark, col_class_eve]

    df_report = df_dark2[cols_report].copy()
    df_report = df_report.sort_values(by=cols_report, ascending=True).reset_index()[cols_report]
    return df_report

def get_classification_report(
    df_dark: pd.DataFrame, path_classifier: str,
    feat_eves: list, x_range_dark: list,
    col_class_dark: str, col_class_eve: str, col_id_dark: str,
    verbose = False
) -> pd.DataFrame:
    ''' 
        TODO: 2023-06-17 - ADD Description
    ''' 

    # Parse input
    X_dark, y_dark, ids = get_input_x_y_ids(
        df_dark=df_dark, x_range_dark=x_range_dark,
        col_id_dark=col_id_dark, col_class_dark=col_class_dark
    )
    n_samples, n_features = X_dark.shape

    
    if (verbose):
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

    # Transform data to classify
    df_dark2, feat_dark2 = get_classification_df(
        X_dark=X_dark, df_dark=df_dark,
        feat_eves=feat_eves,
        col_class_dark=col_class_dark, col_id_dark=col_id_dark
    )

    if (verbose):
        print('-----------------------------------')
        # df_dark2_cp = df_dark2.copy()

        n_match_feats = sum(feat_dark2 == feat_eves)
        print(f'feat_eves: {len(feat_eves)}') 
        print(f'feat_dark2: {len(feat_dark2)}')
        print(f'n_match_feats: {n_match_feats}')

    y_hat = classify(df_dark2=df_dark2, feat_dark2=feat_dark2, path_classifier=path_classifier)
    df_dark2[col_class_eve] = y_hat.copy()

    df_report = build_report(
        # df_dark2=df_dark2, y_hat=y_hat,
        df_dark2=df_dark2,
        col_class_dark=col_class_dark, col_class_eve=col_class_eve, col_id_dark=col_id_dark
    )

    return df_report
    # return df_report, y_hat, df_dark2_cp