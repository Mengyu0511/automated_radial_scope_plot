
def remove_r_from_smi(smi, r_num):
    """ Remove the [*:1] part of a smiles for an r group """

    if smi.endswith(f"[*:{r_num}]"):
        smi = smi[:-5]
    return smi

def get_smiles_and_activity(df, r_groups, activity_col='conversion', to_int=True):
    """
    Returns a dict of r groups, containing the smiles and the activity data for use in the radial scope plot

    r_data will look like - {'R1': {'C': 85, 'Br': 45, ...},
                             'R2": ...}
    """

    r_data = {}

    for num in r_groups:
        smi_col = f"R{num}_smiles"
        r = f"R{num}"

        # empty_dict to hold r group smiles and activity data
        r_data[r] = {}

        # r_df - only rows where the r_group column is not none.
        r_df = df.dropna(subset=[smi_col])

        for index, row in r_df.iterrows():
            # get the smi for this r group
            smi = row[smi_col]
            smi = remove_r_from_smi(smi, num)

            # get the activity for this r group
            activity = row[activity_col]
            if to_int == True:
                activity = int(activity)

            # add smi and activity to dict
            r_data[r][smi] = activity

    return r_data

def combine_enzyme_dict(enzyme_dict_1, enzyme_dict_2):

    combined_data = {k: {} for k in enzyme_dict_1.keys()}

    for r in enzyme_dict_1:
        for smi in enzyme_dict_1[r]:
            combined_data[r][smi] = [enzyme_dict_1[r].get(smi, None),
                                     enzyme_dict_2[r].get(smi, None)]

        for smi in enzyme_dict_2[r]:
            if smi not in combined_data[r]:
                combined_data[r][smi] = [enzyme_dict_1[r].get(smi, None),
                                         enzyme_dict_2[r].get(smi, None)]

    return combined_data