import logging
import logging.config


def setup_logging():
    """Setup basic logging configuration"""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s - %(message)s"
    )
