# -*- coding: utf-8 -*-
#
# test_lfp_recorder.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

"""
Tests for the lfp_recorder
"""

import unittest
import nest


@nest.check_stack
class LFPRecorderTestCase(unittest.TestCase):
    """lfp_recorder tests"""

    def setUp(self):
        nest.ResetKernel()

    def test_SimpleNetwork(self):
        """Test that we record lfp after a spike"""

        lfp = nest.Create('lfp_recorder', 1,
                          {'tau_rise': [1/0.48886579219795934],
                           'tau_decay': [1./0.48856915462005396]})
        s_generator = nest.Create('spike_generator', 1, {'spike_times': [10.]})
        multi = nest.Create('multimeter', 1, {'record_from': ['lfp']})

        nest.Connect(s_generator, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi, lfp)

        nest.Simulate(35.)

        recorded_lfp = nest.GetStatus(multi)[0]['events']['lfp']
        self.assertTrue(len(recorded_lfp) > 0.0)
        # We should not have a LFP signal until after we have had a spike
        self.assertEqual(recorded_lfp[10], 0.0)
        self.assertFalse(recorded_lfp[11] == 0.0)

    def test_SeveralDetectors(self):
        """Test that accumulation of several lfp signals equals lfp signal of
        single detector"""

        # First create several lfp_recorders, spike_generators and multimeters:
        lfp1 = nest.Create('lfp_recorder', 1,
                           {'tau_rise': [1/0.48886579219795934],
                            'tau_decay': [1./0.48856915462005396]})
        s_generator1 = nest.Create('spike_generator', 1,
                                   {'spike_times': [10., 15., 20., 25., 30.,
                                                    35., 40., 45., 50., 55.]})
        multi1 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp2 = nest.Create('lfp_recorder', 1,
                           {'tau_rise': [1/0.48886579219795934],
                            'tau_decay': [1./0.48856915462005396]})
        s_generator2 = nest.Create('spike_generator', 1,
                                   {'spike_times': [12., 17., 22., 27., 32.,
                                                    37., 42., 47., 52., 57.]})
        multi2 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp3 = nest.Create('lfp_recorder', 1,
                           {'tau_rise': [1/0.48886579219795934],
                            'tau_decay': [1./0.48856915462005396]})
        s_generator3 = nest.Create('spike_generator', 1,
                                   {'spike_times': [14., 19., 24., 29., 34.,
                                                    39., 44., 49., 54., 59.]})
        multi3 = nest.Create('multimeter', 1, {'record_from': ['lfp']})

        # Connect:
        nest.Connect(s_generator1, lfp1, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi1, lfp1)
        nest.Connect(s_generator2, lfp2, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi2, lfp2)
        nest.Connect(s_generator3, lfp3, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi3, lfp3)

        # Simulate:
        nest.Simulate(200.)

        # Collect:
        gs1 = nest.GetStatus(multi1)
        gs2 = nest.GetStatus(multi2)
        gs3 = nest.GetStatus(multi3)

        collected_lfp = (gs1[0]['events']['lfp'] +
                         gs2[0]['events']['lfp'] +
                         gs3[0]['events']['lfp'])

        # Now we set up a system with a single lfp_recorder. We still have
        # several spike generators.
        # The LFP signal from this sould equal the accumulation of the LFP
        # signals of the three lfp_recorders above.

        nest.ResetKernel()

        lfp = nest.Create('lfp_recorder', 1,
                          {'tau_rise': [1/0.48886579219795934],
                           'tau_decay': [1./0.48856915462005396]})
        multi = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        s_generator1 = nest.Create('spike_generator', 1,
                                   {'spike_times': [10., 15., 20., 25., 30.,
                                                    35., 40., 45., 50., 55.]})
        s_generator2 = nest.Create('spike_generator', 1,
                                   {'spike_times': [12., 17., 22., 27., 32.,
                                                    37., 42., 47., 52., 57.]})
        s_generator3 = nest.Create('spike_generator', 1,
                                   {'spike_times': [14., 19., 24., 29., 34.,
                                                    39., 44., 49., 54., 59.]})

        nest.Connect(s_generator1, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(s_generator2, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(s_generator3, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi, lfp)

        nest.Simulate(200.)
        gs = nest.GetStatus(multi)

        for count in range(len(collected_lfp)):
            self.assertAlmostEqual(collected_lfp[count],
                                   gs[0]['events']['lfp'][count],
                                   14)

        # Lastly, we test what happens if we have one spike_generator. It
        # should equal the two LFP signals from above.
        nest.ResetKernel()

        lfp = nest.Create('lfp_recorder', 1,
                          {'tau_rise': [1/0.48886579219795934],
                           'tau_decay': [1./0.48856915462005396]})
        multi = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        s_generator = nest.Create('spike_generator', 1,
                                  {'spike_times': [10., 12., 14., 15., 17.,
                                                   19., 20., 22., 24., 25.,
                                                   27., 29., 30., 32., 34.,
                                                   35., 37., 39., 40., 42.,
                                                   44., 45., 47., 49., 50.,
                                                   52., 54., 55., 57., 59.]})

        nest.Connect(s_generator, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi, lfp)

        nest.Simulate(200.)
        gs2 = nest.GetStatus(multi)

        for count in range(len(collected_lfp)):
            self.assertAlmostEqual(collected_lfp[count],
                                   gs2[0]['events']['lfp'][count],
                                   14)
            self.assertAlmostEqual(gs[0]['events']['lfp'][count],
                                   gs2[0]['events']['lfp'][count],
                                   14)

    def test_SeveralDetectorsWithPoissonInput(self):
        """Test lfp_recorder with poisson input"""
        # This is basically the same test as above, but with poisson_generator
        # instead of spike_generators.

        # Create:
        lfp1 = nest.Create('lfp_recorder', 1,
                           {'tau_rise': [1/0.48886579219795934],
                            'tau_decay': [1./0.48856915462005396]})
        p_generator1 = nest.Create('poisson_generator', 1, {'rate': 800.})
        multi1 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp2 = nest.Create('lfp_recorder', 1,
                           {'tau_rise': [1/0.48886579219795934],
                            'tau_decay': [1./0.48856915462005396]})
        p_generator2 = nest.Create('poisson_generator', 1, {'rate': 800.})
        multi2 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp3 = nest.Create('lfp_recorder', 1,
                           {'tau_rise': [1/0.48886579219795934],
                            'tau_decay': [1./0.48856915462005396]})
        p_generator3 = nest.Create('poisson_generator', 1, {'rate': 800.})
        multi3 = nest.Create('multimeter', 1, {'record_from': ['lfp']})

        # Connect
        nest.Connect(p_generator1, lfp1, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi1, lfp1)
        nest.Connect(p_generator2, lfp2, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi2, lfp2)
        nest.Connect(p_generator3, lfp3, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi3, lfp3)

        # Simulate
        nest.Simulate(200.)

        # Collect
        gs1 = nest.GetStatus(multi1)
        gs2 = nest.GetStatus(multi2)
        gs3 = nest.GetStatus(multi3)

        collected_lfp = (gs1[0]['events']['lfp'] +
                         gs2[0]['events']['lfp'] +
                         gs3[0]['events']['lfp'])

        # Single
        nest.ResetKernel()

        lfp = nest.Create('lfp_recorder', 1,
                          {'tau_rise': [1/0.48886579219795934],
                           'tau_decay': [1./0.48856915462005396]})
        multi = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        p_generator1 = nest.Create('poisson_generator', 1, {'rate': 800.})
        p_generator2 = nest.Create('poisson_generator', 1, {'rate': 800.})
        p_generator3 = nest.Create('poisson_generator', 1, {'rate': 800.})

        nest.Connect(p_generator1, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(p_generator2, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(p_generator3, lfp, 'one_to_one',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi, lfp)

        nest.Simulate(200.)
        gs = nest.GetStatus(multi)

        for count in range(len(collected_lfp)):
            self.assertAlmostEqual(collected_lfp[count],
                                   gs[0]['events']['lfp'][count])


def suite():
    suite = unittest.makeSuite(LFPRecorderTestCase, 'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


if __name__ == "__main__":
    run()
