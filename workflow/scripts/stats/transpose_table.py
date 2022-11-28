import pandas as pd
import os


def write_to_html_file(df, title):
    inp = """
        <html>
        <head>
        <style>
            h2 {
                text-align: center;
                font-family: Helvetica, Arial, sans-serif;
            }
            table { 
                margin-left: auto;
                margin-right: auto;
            }
            table, th, td {
                border: 1px solid black;
                border-collapse: collapse;
            }
            th, td {
                padding: 5px;
                text-align: center;
                font-family: Helvetica, Arial, sans-serif;
                font-size: 90%;
            }
            table tbody tr:hover {
                background-color: #dddddd;
            }
            .wide {
                width: 90%; 
            }
        </style>
        </head>
        <body>
        """
    out = """
            </body>
            </html>
            """
    result = inp
    result += "<h2> {} statistics summary </h2>\n".format(str(title))
    result += df.to_html(classes="wide", escape=False, index=False)
    result += out
    return result


df = pd.read_csv(snakemake.input[0], sep="\t")
df["callset"] = df["callset"].apply(lambda r: os.path.basename(r).replace(".tsv", ""))
df = df.set_index("callset")
df = df.fillna(0).T.reset_index()
pd.options.display.float_format = "{:,.1f}".format
df_out = write_to_html_file(df, snakemake.wildcards.sample)
with open(snakemake.output.html, "w") as o:
    o.write(df_out)
