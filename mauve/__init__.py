import os
LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(LOCAL_DIR, '..')))

if not os.environ.get('MAUVE_DIR'):
    try:
        os.environ['MAUVE_DIR'] = os.path.join(LOCAL_DIR, 'src', 'bin')
    except:
        raise ImportError('MAUVE_DIR environmental variable is not set.')

from .buildindex import buildindex, indexutils

__author__ = "Nick Conway, Ben Pruitt"
__all__ = ['buildindex', 'indexutils']