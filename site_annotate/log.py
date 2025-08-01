import os
import logging


# def get_logger(name):

#     # Create a logger object
#     logger = logging.getLogger(name)
#     logger.setLevel(logging.INFO)  # Set the default logging level

#     # Create handlers
#     console_handler = logging.StreamHandler()
#     file_handler = logging.FileHandler("app.log")

#     # Set level for handlers
#     console_handler.setLevel(logging.INFO)
#     file_handler.setLevel(logging.INFO)

#     # Create a logging format
#     formatter = logging.Formatter(
#         "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
#     )
#     console_handler.setFormatter(formatter)
#     file_handler.setFormatter(formatter)

#     # Add handlers to the logger
#     logger.addHandler(console_handler)
#     # logger.addHandler(file_handler)

#     return logger


def get_logger(name: str, log_file: str = "app.log") -> logging.Logger:
    logger = logging.getLogger(name)

    if logger.hasHandlers():
        return logger  # Prevent duplicate handlers if already configured

    # logger.setLevel(logging.INFO)
    logger.setLevel(logging.DEBUG)

    # Console handler
    console_handler = logging.StreamHandler()
    # console_handler.setLevel(logging.INFO)
    console_handler.setLevel(logging.DEBUG)

    # File handler (create log directory if needed)
    log_dir = os.path.dirname(log_file)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setLevel(logging.INFO)

    # Formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Add handlers
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger
