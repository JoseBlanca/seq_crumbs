
import os.path

from crumbs.utils.bin_utils import create_get_binary_path
from bam_crumbs.settings import get_setting


BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                       '..', '..', 'bin'))

get_binary_path = create_get_binary_path(os.path.split(__file__)[0],
                                         get_setting)
