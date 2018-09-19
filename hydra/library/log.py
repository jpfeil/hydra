import logging

def setup_custom_logger(name, filename, level):
    """
    https://stackoverflow.com/questions/7621897/python-logging-module-globally

    :param name:
    :return:
    """
    logging.basicConfig(filename=filename,
                        level=level)

    logging.getLogger().addHandler(logging.StreamHandler())

    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger