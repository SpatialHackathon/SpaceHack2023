#!/bin/bash

# Output config file
CONFIG_FILE="generated_path_config.yml"

# Function to write the YAML structure comment block to a file
cat <<'EOF' > "$CONFIG_FILE"
# The yaml file follows the following structure

# * categories (datasets/methods/metrics)
#     - {name}
#       - env: path/to/conda/env/.yaml
#       - script: path/to/script/.{py|r}
#       - env_additional: (optional)path/to/installation/script/.sh
#       - optargs: path/to/input/parameters/.json

# * config_files (for methods/metrics)
#     - {name} # MUST BE THE SAME AS THE METHOD/METRIC NAME
#       - {config_name}: path/to/config
#       - script: path to the excutation script
#       - env_additional: Only for certain methods, need installation shell script (.sh)
#       - optargs: optional arguments file (for input control/quality control)

# Notice for new addition:
# - name must be the same as the folder name!
# - All identation is 2 spaces!
# - When adding methods/metrics, remember to also add config_files if avaliable! 
# - Comment out configs that you don't want to run.
EOF
echo >> $CONFIG_FILE

# Function to extract the '*' part from a .yml or other .format file in a folder
extract_filename() {
    local path="$1"

    # Find the .yml file in the folder
    local file=$(find "$path" -maxdepth 1 -name "*.$2" -print -quit)

    if [[ -n "$file" ]]; then
        # Extract the basename and remove the .yml extension
        local filename=$(basename "$file")
        local extracted=${filename%.2}
        echo "$extracted"
    else
        echo "none"
    fi
}

# Start the configuration file
echo "datasets:" >> $CONFIG_FILE

# Process datasets
for dataset in ../data/*/; do
    name=$(basename "$dataset")
    dataset_format="${dataset/\.\.\//}"

    yml=$(extract_filename $dataset "yml")
    optargs=$(extract_filename $dataset "json")
    env_additional=$(extract_filename $dataset "sh")

    r_script=$(extract_filename $dataset "r")
    python_script=$(extract_filename $dataset "py")

    script=$r_script

    if [[ $r_script == "none" ]]; then
        script=$python_script
    fi

    echo "  $name:" >> $CONFIG_FILE
    echo "    env: $dataset_format$yml" >> $CONFIG_FILE
    echo "    script: $dataset_format$script" >> $CONFIG_FILE
    echo "    optargs: $dataset_format$optargs" >> $CONFIG_FILE

    if [[ $env_additional != "none" ]]; then
        echo "    env_additional: $dataset_format$env_additional" >> $CONFIG_FILE
    fi
    
    echo >> $CONFIG_FILE
done

echo "neighbors_infos:" >> $CONFIG_FILE
echo "  delaunay_triangulation:" >> $CONFIG_FILE
echo "    script: preprocessing/neighbors/delaunay_triangulation/delaunay_triangulation.py" >> $CONFIG_FILE
echo "    env: preprocessing/neighbors/delaunay_triangulation/delaunay_triangulation.yml" >> $CONFIG_FILE
echo >> $CONFIG_FILE

# Add methods section
echo "methods:" >> $CONFIG_FILE
for method in ../method/*/; do
    name=$(basename "$method")
    method_format="${method/\.\.\//}"

    yml=$(extract_filename $method "yml")
    optargs=$(extract_filename $method "json")
    env_additional=$(extract_filename $method "sh")

    r_script=$(extract_filename $method "r")
    python_script=$(extract_filename $method "py")

    script=$r_script

    if [[ $r_script == "none" ]]; then
        script=$python_script
    fi

    echo "  $name:" >> $CONFIG_FILE
    echo "    env: $method_format$yml" >> $CONFIG_FILE
    echo "    script: $method_format$script" >> $CONFIG_FILE
    echo "    optargs: $method_format$optargs" >> $CONFIG_FILE

    if [[ $env_additional != "none" ]]; then
        echo "    env_additional: $method_format$env_additional" >> $CONFIG_FILE
    fi
    
    echo >> $CONFIG_FILE
done

# Add metrics section
echo "metrics:" >> $CONFIG_FILE
for metric in ../metric/*/; do
    name=$(basename "$metric")
    metric_format="${metric/\.\.\//}"

    yml=$(extract_filename $metric "yml")
    optargs=$(extract_filename $metric "json")
    env_additional=$(extract_filename $metric "sh")

    r_script=$(extract_filename $metric "r")
    python_script=$(extract_filename $metric "py")

    script=$r_script

    if [[ $r_script == "none" ]]; then
        script=$python_script
    fi

    echo "  $name:" >> $CONFIG_FILE
    echo "    env: $metric_format$yml" >> $CONFIG_FILE
    echo "    script: $metric_format$script" >> $CONFIG_FILE
    echo "    optargs: $metric_format$optargs" >> $CONFIG_FILE

    if [[ $env_additional != "none" ]]; then
        echo "    env_additional: $metric_format$env_additional" >> $CONFIG_FILE
    fi
    
    echo >> $CONFIG_FILE
done

# Add config_files section
echo "config_files:" >> $CONFIG_FILE
for method in ../method/*/; do
    name=$(basename "$method")

    echo "  $name:" >> $CONFIG_FILE

    # Check if config_defaul is present
    if [[ ! -f "$method/config/config_default.json" ]]; then
        echo "Error: Required file '$method/config/config_default.json' does not exist in the folder." >&2
        exit 1
    fi

    for config in "$method"/config/*.json; do
        configname=$(basename "$config")
        configname="${configname/\.json/}"
        echo "    $configname: config/$configname.json" >> $CONFIG_FILE
    done
done

echo >> $CONFIG_FILE

# Notify completion
echo "Configuration file generated at $CONFIG_FILE."