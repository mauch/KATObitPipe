# Module to assist in running AIPS/ParselTongue on machines without
# an AIPS installation.
#
# Stephen Bourke, JIVE
# Version: 20090204

"""
This module assists in running AIPS/ParselTongue on machines without
an AIPS installation.
"""

import os, shutil, string, sys, time

aips_server = 'ftp.aoc.nrao.edu'

default_year = time.gmtime().tm_year - 1
default_version = '31DEC' + str(default_year)[-2:]

# Minimum files required:
intel_libs = ['libimf.so', 'libsvml.so']
mac_libs = ['libimf.dylib', 'libirc.dylib', 'libsvml.dylib']
popsdat_files = ['POPSDAT.HLP']
binary_files = ['FILAIP.EXE']

def arch():
    """Return the AIPS architecture for the current machine"""
    uname = os.uname()
    if uname[0] == 'Linux' and ((uname[4][0] == 'i' and uname[4][-2:] == '86')):
        return 'LINUX'
    elif uname[0] == 'Linux' and (uname[4] == 'x86_64'):
        return 'LNX64'
    elif uname[0] == 'Darwin':
        return 'MACINT'
    else:
        raise NotImplementedError('Unknown Architecture')

def init_environ(path=None, version=default_version):
    """Set required environment variables"""
    if not path:
        path = os.getcwd()
    os.environ['AIPS_ROOT'] = path
    os.environ['AIPS_VERSION'] = os.environ['AIPS_ROOT'] + '/' + version
    os.environ['ARCH'] = arch()
    os.environ['LOAD'] = os.environ['AIPS_VERSION'] + '/' + os.environ['ARCH'] + '/LOAD'
    os.environ['NEW'] = os.environ['AIPS_VERSION']
    os.environ['VERSION'] = 'NEW'
    os.environ['PLOTTER'] = '/tmp'
    # Add the intel libs to [DY]LD_LIBRARY_PATH
    if (os.environ['ARCH'] == 'LINUX') or (os.environ['ARCH'] == 'LNX64'):
        lib_env = 'LD_LIBRARY_PATH'
    elif os.environ['ARCH'] == 'MACINT':
        lib_env = 'DYLD_LIBRARY_PATH'
    lib_dir = '%s/%s/LIBR/INTELCMP' % (os.environ['AIPS_VERSION'], os.environ['ARCH'])
    if lib_env in os.environ:
        os.environ[lib_env] += ':' + lib_dir
    else:
        os.environ[lib_env] = lib_dir

def create_path_list(path, file_list):
    url_list = []
    for filename in file_list:
        url = '%s/%s' % (path, filename)
        url_list.append(url)
    return url_list

def lib_urls():
    base = version() + '/' + os.environ['ARCH'] + '/LIBR/INTELCMP'
    if (os.environ['ARCH'] == 'LINUX') or (os.environ['ARCH'] == 'LNX64'):
        lib_files = intel_libs
    elif os.environ['ARCH'] == 'MACINT':
        lib_files = mac_libs
    return create_path_list(base, lib_files)

def popsdat_urls():
    base = version() + '/HELP'
    return create_path_list(base, popsdat_files)

def binary_urls():
    exe_path = version() + '/' + os.environ['ARCH'] + '/LOAD'
    return create_path_list(exe_path, binary_files)

def rsync(server, path_list, force=False):
    args = ['rsync', '--compress', '--relative', '--no-motd', '--progress']
    if not force:
        args.append('--ignore-existing')
    args.append(server + '::' + ' '.join(path_list))
    args.append(os.environ['AIPS_VERSION'])    # Download destination
    if os.spawnvp(os.P_WAIT, 'rsync', args) != 0:
        # TODO: Better failure reporting
        print('rsync failure', file=sys.stderr)
        sys.exit(-1)

def make_da00(da00_path=None, force=False):
    """Create DA00 directory."""
    if not da00_path:
        da00_path = os.environ['AIPS_VERSION'] + '/DA00'
    if not os.path.exists(da00_path) or force:
        template_path = os.environ['AIPS_VERSION'] + '/' + os.environ['ARCH'] + '/TEMPLATE'
        shutil.copytree(template_path, da00_path)
    os.environ['DA00'] = da00_path
    os.environ['NET0'] = da00_path
    # Added by TM 25/3/2012 Set up NETSP file in $NET0 to remove error messages.
    open(da00_path+'/NETSP', 'w').close()

def make_disk(disk_path=None):
    """Make the directory and put the 'SPACE' file in it."""
    if not disk_path:
        disk_path = os.environ['AIPS_VERSION'] + '/DATA'
    space_file = disk_path + '/SPACE'
    if not os.path.exists(space_file):
        if not os.path.isdir(disk_path):
            os.makedirs(disk_path)
        open(space_file, 'w').close()

    disk_list_add(disk_path)

def ehex(num, width=0, pad_char='0'):
    """Convert a number to base 36."""
    chars = string.digits + string.ascii_uppercase
    base = len(chars)
    hex_str = ''
    while num > 0:
        hex_str = chars[num % base] + hex_str
        num //= base
    pad_len = width - len(hex_str)
    if pad_len > 0:
        hex_str = str(pad_char) * pad_len + hex_str
    return hex_str

def disk_list_add(*args):
    """Add DAxx entries for the directories supplied and update NVOL."""
    try:
        nvol = os.environ['NVOL']
    except KeyError:
        nvol = 0
    for i in range(len(args)):
        os.environ['DA%s' % ehex(nvol+1+i, 2, 0)] = args[i]
    os.environ['NVOL'] = str(nvol + len(args))

def get_aips(basedir=None, version=default_version, force=False):
    """Get all the required files from NRAO and setup the environment"""
    if not basedir:
        basedir = os.getcwd()
    init_environ(basedir, version=version)
    urls = []
    for url in lib_urls() + popsdat_urls() + binary_urls():
        if not os.path.exists(os.environ['AIPS_ROOT'] +'/' + url):
            urls = urls + [url]
    if len(urls)>0:
        rsync(aips_server, urls, force=force)
    filaip(force=force)

def filaip(force=False,data_dir=None):
    mem_dir = os.environ['AIPS_VERSION'] + '/' + os.environ['ARCH'] + '/MEMORY'
    template_dir = os.environ['AIPS_VERSION'] + '/' + os.environ['ARCH'] + '/TEMPLATE'
    if data_dir==None: data_dir = os.environ['AIPS_VERSION'] + '/DATA'
    if force:
        shutil.rmtree(mem_dir, ignore_errors=True)
        shutil.rmtree(template_dir, ignore_errors=True)
        msfile = data_dir + '/MSD001000.001;'
        if os.path.exists(msfile):
            os.remove(msfile)
    if os.path.exists(mem_dir) or os.path.exists(template_dir):
        run_filaip = False
    else:
        run_filaip = True
    for tdir in [mem_dir, template_dir, data_dir]:
        if not os.path.exists(tdir):
            os.makedirs(tdir)
    # Download DA00 to TEMPLATE. make_da00 will copy this to the actual DA00
    os.environ['DA00'] = template_dir
    os.environ['NET0'] = template_dir
    os.environ['DA01'] = data_dir
    os.environ['NVOL'] = '1'
    os.environ['NEWMEM'] = mem_dir
    if run_filaip:
        os.system('echo 8 2 | %s/%s/LOAD/FILAIP.EXE' % (os.environ['AIPS_VERSION'], os.environ['ARCH']))
    for var in ['DA00', 'NET0', 'NVOL']:
        del os.environ[var]

def setup_all(basedir=None, version=default_version, force=False):
    """Get required files and make DA00 and DISK areas"""
    get_aips(basedir=basedir, version=version, force=force)
    make_da00(force=force)
    make_disk()

def version():
    """Return the AIPS version"""
    return os.environ['AIPS_VERSION'].split('/')[-1]

def get_task(*args, **kwargs):
    """Get the Task and Help file for the specified task"""
    try:
        force = kwargs['force']
    except KeyError:
        force = False
    exe_path = version() + '/' + os.environ['ARCH'] + '/LOAD'
    help_path = version() + '/HELP'
    urls = []
    # TM 27/3/13
    # Modified to not go to rsync if the file already exists to speed things up.
    for taskname in args:
        exename = exe_path +'/'+ taskname.upper() + '.EXE'
        hlpname = help_path +'/'+ taskname.upper() + '.HLP'
        if not (os.path.exists(os.environ['AIPS_ROOT'] + '/' + exename) or os.path.exists(os.environ['AIPS_ROOT'] + '/' + hlpname)):
            urls += create_path_list(exe_path, [taskname.upper() + '.EXE'])
            urls += create_path_list(help_path, [taskname.upper() + '.HLP'])
    if len(urls)>0:
        rsync(aips_server, urls, force=force)

quick_start = setup_all
