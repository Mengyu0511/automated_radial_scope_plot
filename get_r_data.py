

def get_smiles_and_activity(df, r_groups, activity_col='conversion', to_int=True):
    """
    Returns two dicts containing lists of the smiles and the activity data for use in the radial scope plot
    """

    r_activity = {}
    r_smiles = {}

    for num in r_groups:
        column = f"R{num}_smiles"
        r = f"R{num}"

        # r_df - only rows where the r_group column is not none.
        r_df = df.dropna(subset=[column])

        # get activity for those rows in r_df
        activity_list = r_df[activity_col].tolist()
        if to_int == True:
            activity_list = [int(c) for c in activity_list]
        r_activity[r] = activity_list

        # get smiles for those rows
        r_list = r_df[column].dropna().tolist()
        for i, sml in enumerate(r_list):
            if sml.endswith(f"[*:{num}]"):
                r_list[i] = r_list[i][:-5]

        r_smiles[r] = r_list

    return r_smiles, r_activity

