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
import math
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

    def test_kernel_convolution(self):
        """Test lfp with realistic data"""

        # To test that the lfp we get from nestkernel is realistic, we test it
        # against a convolution between the firing rate and a kernel on python
        # level, where the kernel is known to give realistic values.
        # The kernel will only fit with the tau_rise, tau_decay and normalizer
        # given below, as those values are fitted to this kernel.
        kernel = [1.10167120e-16, -3.21767675e-16, 8.50309087e-16,
                  -2.10825843e-15, 4.98932743e-15, -1.13729372e-14,
                  2.50921741e-14, -5.37078004e-14, 1.11562884e-13,
                  -2.24559409e-13, 4.36342190e-13, -8.12519563e-13,
                  1.43109033e-12, -2.32320376e-12, 3.23745828e-12,
                  -3.04445770e-12, -9.35353620e-13, 2.13181578e-11,
                  -1.10221233e-10, 1.13624809e-10, 2.63285290e-11,
                  9.84213201e-09, 1.30181595e-08, 1.13629411e-08,
                  8.39004294e-09, 5.30052880e-09, 3.57433098e-09,
                  2.20725487e-09, 1.44674216e-09, 8.93289079e-10,
                  5.70139073e-10, 3.51097967e-10, 2.20241171e-10,
                  1.35573201e-10, 8.42163899e-11, 5.18446572e-11,
                  3.20279475e-11, 1.97108533e-11, 1.21370068e-11,
                  7.46330397e-12, 4.58550766e-12]

        # Create, Connect, Simulate
        tau_rise = [1 / 0.4137949562815958]
        tau_decay = [1 / 0.9631522808240071]
        normalizer = [2.6512177702898558e-08]

        lfp = nest.Create('lfp_detector', 1, {'tau_rise': tau_rise,
                                              'tau_decay': tau_decay,
                                              'normalizer': normalizer})
        parrot = nest.Create('parrot_neuron', 50)
        multi = nest.Create('multimeter', 1, {'record_from': ['lfp']})
        spikes = nest.Create('spike_detector')
        p_generator = nest.Create('poisson_generator', 1, {'rate': 400.})

        nest.Connect(p_generator, parrot)
        nest.Connect(parrot, lfp, 'all_to_all',
                     {'model': 'static_synapse', 'receptor_type': 1})
        nest.Connect(multi, lfp)
        nest.Connect(parrot, spikes)

        simtime = 600.
        nest.Simulate(simtime)

        lfp_nest = nest.GetStatus(multi)[0]['events']['lfp']

        # Get spiking rate
        spike_times = nest.GetStatus(spikes)[0]['events']['times']
        rate = [0] * int(simtime-1)
        for t in spike_times:
            t_rounded = math.floor(t)
            rate[t_rounded-1] += 1

        # Convolve kernel with the spiking rate to get realistic values
        lfp_convolution = np.convolve(rate, kernel, 'same')

        for count in range(len(lfp_nest)):
            # lfp is given in mV, and should not be bigger than 1 mV.
            self.assertTrue(abs(lfp_nest[count]) < 1.)
            self.assertAlmostEqual(lfp_nest[count],
                                   lfp_convolution[count], 6)

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

        normalizer = [0.00015807504272622632,
                      1.1871740309906909e-05,
                      -6.241511286160071e-06,
                      4.84509284785128e-07,
                      -9.825292432345066e-06,
                      2.10407553698354e-06,
                      -1.6348916609249954e-05,
                      1.8102033694846387e-06,
                      -0.002474980238314617]

        # First we get individual LFP contribution
        borders = [n for i in range(5, 8) for n in [i, i]]
        pg = nest.Create('poisson_generator', 1, {'rate': 80000.})
        lfp = [nest.Create('lfp_detector', 1,
                           {'tau_rise': tau_rise,
                            'tau_decay': tau_decay,
                            'normalizer': normalizer,
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
            self.assertFalse(np.sum(individual) == 0.)

        # Manual sum of LFP contributions
        sum_of_recorded_lfp = np.sum(recorded_lfp, axis=0)

        # Then we get the summed LFP contribution
        nest.ResetKernel()
        borders = [n for i in range(3, 6) for n in [i, i]]
        pg = nest.Create('poisson_generator', 1, {'rate': 80000.})
        lfp = nest.Create('lfp_detector', 1,
                          {'tau_rise': tau_rise,
                           'tau_decay': tau_decay,
                           'normalizer': normalizer,
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
        self.assertFalse(np.sum(recorded_lfp_sum) == 0.)

        # LFP values should be (almost) equal
        np.testing.assert_array_almost_equal(sum_of_recorded_lfp,
                                             recorded_lfp_sum)


def suite():
    suite = unittest.makeSuite(LFPDetectorTestCase, 'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


if __name__ == "__main__":
    run()
