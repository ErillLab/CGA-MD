import logging

"""
Documentation:

This script configures the logging module to log messages to a file named 'log.log' in append mode.
Log messages will include the timestamp, log level, and the message itself.

Usage:
    Ensure this script is imported or executed before any logging calls are made in the application.
"""

# Configure logging settings
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                    filename="log.log",
                    filemode="w")

# Create a logger
logger = logging.getLogger(__name__)

# Add a handler for error messages
error_handler = logging.FileHandler("log.log")
error_handler.setLevel(logging.ERROR)
error_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
error_handler.setFormatter(error_formatter)
logger.addHandler(error_handler)