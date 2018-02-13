# -*- coding: utf-8 -*-
#
# test_lfp_detector.py
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
Tests for the lfp_detector
"""

import unittest
import nest
import numpy as np


@nest.check_stack
class LFPDetectorTestCase(unittest.TestCase):
    """lfp_detector tests"""

    def setUp(self):
        nest.ResetKernel()

    def test_SimpleNetwork(self):
        """Test that we record lfp after a spike"""

        lfp = nest.Create('lfp_detector', 1,
                          {'tau_rise': [1 / 0.48886579219795934],
                           'tau_decay': [1. / 0.48856915462005396],
                           'normalizer': [0.00015807504272622632]})
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

        # First create several lfp_detectors, spike_generators and multimeters:
        lfp1 = nest.Create('lfp_detector', 1,
                           {'tau_rise': [1 / 0.48886579219795934],
                            'tau_decay': [1. / 0.48856915462005396],
                            'normalizer': [0.00015807504272622632]})
        s_generator1 = nest.Create('spike_generator', 1,
                                   {'spike_times': [10., 15., 20., 25., 30.,
                                                    35., 40., 45., 50., 55.]})
        multi1 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp2 = nest.Create('lfp_detector', 1,
                           {'tau_rise': [1 / 0.48886579219795934],
                            'tau_decay': [1. / 0.48856915462005396],
                            'normalizer': [0.00015807504272622632]})
        s_generator2 = nest.Create('spike_generator', 1,
                                   {'spike_times': [12., 17., 22., 27., 32.,
                                                    37., 42., 47., 52., 57.]})
        multi2 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp3 = nest.Create('lfp_detector', 1,
                           {'tau_rise': [1 / 0.48886579219795934],
                            'tau_decay': [1. / 0.48856915462005396],
                            'normalizer': [0.00015807504272622632]})
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

        # Now we set up a system with a single lfp_detector. We still have
        # several spike generators.
        # The LFP signal from this sould equal the accumulation of the LFP
        # signals of the three lfp_recorders above.

        nest.ResetKernel()

        lfp = nest.Create('lfp_detector', 1,
                          {'tau_rise': [1 / 0.48886579219795934],
                           'tau_decay': [1. / 0.48856915462005396],
                           'normalizer': [0.00015807504272622632]})
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

        lfp = nest.Create('lfp_detector', 1,
                          {'tau_rise': [1 / 0.48886579219795934],
                           'tau_decay': [1. / 0.48856915462005396],
                           'normalizer': [0.00015807504272622632]})
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
        """Test lfp_detector with poisson input"""
        # This is basically the same test as above, but with poisson_generator
        # instead of spike_generators.

        # Create:
        lfp1 = nest.Create('lfp_detector', 1,
                           {'tau_rise': [1 / 0.48886579219795934],
                            'tau_decay': [1. / 0.48856915462005396],
                            'normalizer': [0.00015807504272622632]})
        p_generator1 = nest.Create('poisson_generator', 1, {'rate': 800.})
        multi1 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp2 = nest.Create('lfp_detector', 1,
                           {'tau_rise': [1 / 0.48886579219795934],
                            'tau_decay': [1. / 0.48856915462005396],
                            'normalizer': [0.00015807504272622632]})
        p_generator2 = nest.Create('poisson_generator', 1, {'rate': 800.})
        multi2 = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        lfp3 = nest.Create('lfp_detector', 1,
                           {'tau_rise': [1 / 0.48886579219795934],
                            'tau_decay': [1. / 0.48856915462005396],
                            'normalizer': [0.00015807504272622632]})
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

        lfp = nest.Create('lfp_detector', 1,
                          {'tau_rise': [1 / 0.48886579219795934],
                           'tau_decay': [1. / 0.48856915462005396],
                           'normalizer': [0.00015807504272622632]})
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

    def test_three_population_network(self):
        """LFP is sum invariant"""
        N = 3  # Number of neurons
        T = 100.  # Simulation time
        tau_rise = [1 / a for a in [0.48886579219795934,
                                    0.9035914068978077,
                                    0.1930894552692974,
                                    0.9093497442449103,
                                    0.41639118548249265,
                                    0.9034399617801607,
                                    0.2617577871858636,
                                    0.9059096122950062,
                                    0.710567063977332]]
        tau_decay = [1 / b for b in [0.48856915462005396,
                                     1.3469868506112692,
                                     0.32806051756797333,
                                     1.3498529711342273,
                                     0.41132039426568,
                                     1.3468658588833944,
                                     0.26070589459724114,
                                     1.3481709833572588,
                                     0.7107727064429461]]

        # First we get individual LFP contribution
        borders = [n for i in range(5, 8) for n in [i, i]]
        pg = nest.Create('poisson_generator', 1, {'rate': 80000.})
        lfp = [nest.Create('lfp_detector', 1,
                           {'tau_rise': tau_rise,
                            'tau_decay': tau_decay,
                            'borders': borders}) for i in range(N)]

        neurons = [nest.Create('iaf_psc_alpha') for i in range(N)]
        mm = [nest.Create('multimeter', 1, {'record_from': ['lfp']})
              for i in range(N)]

        for i in range(N):
            nest.Connect(pg, neurons[i])
            nest.Connect(neurons[i], lfp[i], 'one_to_one',
                         {'model': 'static_synapse', 'receptor_type': 1})
            nest.Connect(mm[i], lfp[i])
            for j in range(N):
                nest.Connect(neurons[i], neurons[j])

        nest.Simulate(T)

        recorded_lfp = [nest.GetStatus(multi)[0]['events']['lfp']
                        for multi in mm]
        # Check that every LFP detector recorded something
        for individual in recorded_lfp:
            self.assertTrue(np.sum(individual) > 0.)

        # Manual sum of LFP contributions
        sum_of_recorded_lfp = np.sum(recorded_lfp, axis=0)

        # Then we get the summed LFP contribution
        nest.ResetKernel()
        borders = [n for i in range(3, 6) for n in [i, i]]
        pg = nest.Create('poisson_generator', 1, {'rate': 80000.})
        lfp = nest.Create('lfp_detector', 1,
                          {'tau_rise': tau_rise,
                           'tau_decay': tau_decay,
                           'borders': borders})
        neurons = [nest.Create('iaf_psc_alpha') for i in range(N)]
        mm = nest.Create('multimeter', 1, {'record_from': ['lfp']})

        for i in range(N):
            nest.Connect(pg, neurons[i])
            nest.Connect(neurons[i], lfp, 'one_to_one',
                         {'model': 'static_synapse', 'receptor_type': 1})
            for j in range(N):
                nest.Connect(neurons[i], neurons[j])
        nest.Connect(mm, lfp)

        nest.Simulate(T)

        recorded_lfp_sum = nest.GetStatus(mm)[0]['events']['lfp']
        # Check that LFP detector recorded something
        self.assertTrue(np.sum(recorded_lfp_sum) > 0.)

        # LFP values should be (almost) equal
        np.testing.assert_array_almost_equal(sum_of_recorded_lfp,
                                             recorded_lfp_sum)


def suite():
    suite = unittest.makeSuite(LFPRecorderTestCase, 'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


if __name__ == "__main__":
    run()
