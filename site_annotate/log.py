import logging


def get_logger(name):

    # Create a logger object
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)  # Set the default logging level

    # Create handlers
    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler("app.log")

    # Set level for handlers
    console_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.INFO)

    # Create a logging format
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(console_handler)
    # logger.addHandler(file_handler)

    return logger
