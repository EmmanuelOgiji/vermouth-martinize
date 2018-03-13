# Copyright 2018 University of Groningen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.import pytest

import pytest
import itertools
import os
import vermouth
import vermouth.dssp.dssp as dssp
from vermouth.pdb.pdb import read_pdb
from vermouth.tests.datafiles import (
    PDB_PROTEIN,
    PDB_NOT_PROTEIN,
    PDB_PARTIALLY_PROTEIN,
    DSSP_OUTPUT,
)

DSSP_EXECUTABLE = os.environ.get('MART2_TEST_DSSP', 'dssp')
SECSTRUCT_1BTA = list('CEEEEETTTCCSHHHHHHHHHHHHTCCTTCCCSHHHHHHHHTTT'
                      'SCSSEEEEEESTTHHHHTTTSSHHHHHHHHHHHHHTTCCEEEEEC')


def test_read_dssp2():
    with open(DSSP_OUTPUT) as infile:
        secondary_structure = dssp.read_dssp2(infile)
    assert secondary_structure == SECSTRUCT_1BTA


@pytest.mark.parametrize('savefile', [True, False])
def test_run_dssp(savefile, tmpdir):
    # The test runs twice, once with the savefile set to True so we test with
    # savinf the DSSP output to file, and once with savefile set t False so we
    # do not generate the file. The "savefile" argument is set by
    # pytest.mark.parametrize.
    # The "tmpdir" argument is set by pytest and is the path to a temporary
    # directory that exists only for one iteration of the test.
    if savefile:
        path = tmpdir.join('dssp_output')
    else:
        path = None
    system = vermouth.System()
    system.add_molecule(read_pdb(PDB_PROTEIN))
    secondary_structure = dssp.run_dssp(system,
                                        executable=DSSP_EXECUTABLE,
                                        savefile=path)

    # Make sure we produced the expected sequence of secondary structures
    assert secondary_structure == SECSTRUCT_1BTA

    # If we test with savefile, then we need to make sure the file is created
    # and its content corresponds to the reference (excluding the first lines
    # that are variable or contain non-essencial data read from the PDB file).
    # If we test without savefile, then we need to make sure the file is not
    # created.
    if savefile:
        assert os.path.exists(path)
        with open(path) as genfile, open(DSSP_OUTPUT) as reffile:
            gen = '\n'.join(genfile.readlines()[6:])
            ref = '\n'.join(reffile.readlines()[6:])
            assert gen == ref
    else:
        # Is the directory empty?
        assert not os.listdir(tmpdir)
