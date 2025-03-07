from cerberus import Validator
from logger import logger
import json

def read_config(pathConfig):
    """
    Read configuration settings from a JSON file.

    Parameters:
        pathConfig (str): the path to the JSON config file.

    Returns:
        dict: a dictionary containing the configuration settings.

    Raises:
        FileNotFoundError: if the specified file path does not exist.
        Exception: if the configuration file is invalid according to the defined schema.

    Note:
        This function internally calls `validate_config` to ensure the configuration file adheres to a specific schema.
    """
    
    # Read JSON file
    print("Reading config file...")
    logger.info("Reading config file...")
    with open(pathConfig, "r") as file:
        config = json.load(file)
    # Validating JSON
    print("|--> Validating...")
    logger.info("|--> Validating...")
    if validate_config(config):
        return config

def validate_config(config):
    """
    Validate a config dictionary against a predefined schema.

    Parameters:
        config (dict): the configuration dictionary to validate.

    Returns:
        bool: True if the config is valid according to the schema, False otherwise.

    Raises:
        Exception: if the configuration does not adhere to the predefined schema.

    Note:
        This function uses a predefined JSON schema to validate the configuration dictionary.
    """
    
    # JSON schema
    schema = {
        "entrez_parameters": {
            "type": "dict",
            "schema": {
                "request_limit": {"type": "integer", "min": 1},
                "sleep_time": {"type": "number", "min": 0},
                "email": {"type": "string", "regex": r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$"},
                "api_key": {"type": "string", "nullable": True}
            }
        },
        "blast_parameters": {
            "type": "dict",
            "schema": {
                "db": {"type": "string", "nullable": True},
                "tax_IDs": {"type": "list", "nullable": True},
                "e_value": {"type": "number", "nullable": True},
                "query_coverage": {"type": "number", "nullable": True},
                "max_hits": {"type": "integer", "min": 1, "nullable": True}
            }
        },
        "sequences_parameters": {
            "type": "dict",
            "schema": {
                "max_intergenic_size": {"type": "number", "nullable": False},
                "min_intergenic_size": {"type": "number", "nullable": False},
                "upstream_size_region": {"type": "number", "nullable": False},
                "downstream_size_region": {"type": "number", "nullable": False},
                "max_sequence_length": {"type": "integer", "nullable": True},
                "crawl_back": {"type": "boolean", "nullable": False},
                "fuse_intergenic": {"type": "boolean", "nullable": False}
            }
        },
        "max_identity": {"type": "float", "nullable": False, "min": 0.0, "max": 1.0},
        "meme_parameters": {
            "type": "dict",
            "schema": {
                "mod": {"type": "string", "nullable": False},
                "nmotifs": {"type": "integer", "nullable": False},
                "minw": {"type": "integer", "nullable": False, "min": 0},
                "maxw": {"type": "integer", "nullable": False},
                "revcomp": {"type": "boolean", "nullable": False},
                "pal": {"type": "boolean", "nullable": False}
            }
        },
        "run_only_meme": {"type": "boolean", "nullable": False},
        "output_parameters": {
            "type": "dict",
            "schema": {
                "folder_name": {"type": "string", "nullable": False},
                "folder_rerun_meme": {"type": "string", "nullable": True}
            }
        },
        "input_records": {"type": "list", "minlength": 1, "nullable": False}
    }

    # Validate JSON file
    validator = Validator(schema)
    is_valid = validator.validate(config)

    # Show errors
    if not is_valid:
        for field, error in validator.errors.items():
            raise Exception(f"Error in {field}: {error}")
    else:
        logger.info("|--> JSON file is valid.")
        print("|--> JSON file is valid.")
        return True