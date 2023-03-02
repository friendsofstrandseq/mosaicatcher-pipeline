rule: 

checkpoint filter_input_subclonality:
    container: None
    input: "{folder}/{sample}/scNOVA_input_user/input_subclonality.txt"
    output: "{folder}/{sample}/scNOVA_input_user/input_subclonality_{clone}.txt"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        [df.loc[df["Subclonality"] == clone].to_csv(output[0]) for clone in df.Subclonality.unique().tolist()]
