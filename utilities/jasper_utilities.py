"""
    jasper_utilities.py (BETA)
    Miscellaneous utilities and functions used by the Jasper package.
    This software is a "Camilo Valdes Work" under the terms of the United States Copyright Act. Please cite
    the author(s) in any work or product based on this material. A prototype version of this script was
    created as part of my doctoral work at the BioRG lab at Florida International University, and it was
    expanded at the Quantitative Life Sciences Initiative at the University of Nebraska-Lincoln.
"""

import sys

__author__ = 'Camilo Valdes'
__VERSION__ = 'B1'
__BUILD_NUMBER__ = '2023021075'


def get_jasper_version():
    """
    Returns the current version number of the Jasper framework.
    """
    return __VERSION__


def get_jasper_build_number():
    """
    Returns the current build number of the Jasper framework.
    """
    return __BUILD_NUMBER__


def print_jasper_header():
    """"
        Prints a nice-looking header for displaying in a Terminal.
        Jasper header made with MonoDraw on macOS.
        https://monodraw.helftone.com
        :)
    """

    sys.stdout.flush()

    print(" ───────────────────────────────────────────────────────────────────────────────")
    print("   _____ _______ _______  _____  _______  ______        Version " + get_jasper_version() +
          "." + get_jasper_build_number() + "")
    print("     |   |_____| |______ |_____] |______ |_____/                Jasper")
    print("   __|   |     | ______| |       |______ |    \_         ______(BETA)_______")
    print("                                                                              ")
    print("   https://github.com/camilo-v/jasper                               ")
    print(" ───────────────────────────────────────────────────────────────────────────────")


