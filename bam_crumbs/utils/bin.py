
import os.path
import platform

from crumbs.exceptions import MissingBinaryError
from bam_crumbs.settings import get_setting


BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                       '..', '..', 'bin'))


def get_binary_path(binary_name):
    '''It return the path to the proper binary. It looks on platform and
    architecture to decide it.

    Fails if there is not binary for that architecture
    '''
    if get_setting('USE_EXTERNAL_BIN_PREFIX'):
        binary_name = get_setting('EXTERNAL_BIN_PREFIX') + binary_name

    if not get_setting('ADD_PATH_TO_EXT_BIN'):
        # I have to check if the bynary is on my current directory.
        # If it is there use it, else assumes that it is on the path
        if os.path.exists(os.path.join(os.getcwd(), binary_name)):
            return os.path.join(os.getcwd(), binary_name)
        return binary_name

    system = platform.system().lower()
    if system == 'windows':
        binary_name += '.exe'
    arch = platform.architecture()[0]

    join = os.path.join

    module_path = os.path.split(__file__)[0]
    third_party_path = join(module_path, '..', 'third_party', 'bin')
    third_party_path = os.path.abspath(third_party_path)

    if not os.path.exists(third_party_path):
        msg = 'Third party bin directory not found, please fix me.'
        raise MissingBinaryError(msg)

    binary_path = os.path.abspath(join(third_party_path, system, arch,
                                       binary_name))

    if os.path.exists(binary_path):
        return binary_path
    elif arch == '64bit':
        arch = '32bit'
        binary_path = os.path.abspath(join(third_party_path, system, arch,
                                           binary_name))
        if os.path.exists(binary_path):
            return binary_path

    # At this point there is not available binary for the working platform
    msg = '{} not available for this platform: {}'.format(binary_name, system)
    raise MissingBinaryError(msg)
