import logging

__all__ = ['log']

# Initialize a non-root logger that other files can use
log = logging.getLogger('specutils')
