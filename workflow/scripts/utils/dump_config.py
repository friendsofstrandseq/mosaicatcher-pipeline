import yaml


def update_config(input_file, output_file):
    # Load the existing config file
    with open(input_file, "r") as file:
        flat_file_config = yaml.safe_load(file)

    # Update the config with Snakemake parameters
    for key, value in snakemake.config.items():
        flat_file_config[key] = value

    # Save the updated config to the output file
    with open(output_file, "w") as file:
        yaml.dump(flat_file_config, file)


if __name__ == "__main__":
    input_config = snakemake.input[0]
    output_config = snakemake.output[0]

    update_config(input_config, output_config)
