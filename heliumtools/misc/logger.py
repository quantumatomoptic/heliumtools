import colorlog
import logging


def configure_logger():
    logger = colorlog.getLogger()
    logger.setLevel(logging.INFO)
    

    formatter = colorlog.ColoredFormatter(
        "%(log_color)s%(levelname)s:%(message)s%(reset)s",
        log_colors={
            "DEBUG": "cyan",
            "INFO": "green",
            "WARNING": "yellow",
            "ERROR": "red",
            "CRITICAL": "red,bg_white",
        },
        secondary_log_colors={},
        style="%",
    )

    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)


def getLogger(name, level="info"):
    logger = colorlog.getLogger()
    if level.lower() in "information":
        logger.setLevel(logging.INFO)
    elif level.lower() in "debug no bug":
        logger.setLevel(logging.DEBUG)
    elif level.lower() in "warning error":
        logger.setLevel(logging.WARNING)
    else:
        msg = f"{level} is not a logging level. Setting to info."
        log.error(msg)
        logger.setLevel(logging.DEBUG)
    return logging.getLogger(name)


# Appel de la fonction pour configurer le logger
configure_logger()


if __name__ == "__main__":
    log = getLogger(__name__)
    log.debug("Ceci est un message de débogage.")
    log.info("Ceci est un message informatif.")
    log.warning("Ceci est un avertissement.")
    log.error("Ceci est une erreur.")
    log.critical("Ceci est une erreur critique.")
