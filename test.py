import os
import shutil
import unittest
from mapcounter import MapCounter
from readcounter import ReadCounter
from gccounter import GCCounter


class TestReadCounter(unittest.TestCase):

    def setUp(self):

        ref = 'ref'
        self.tempdir = 'tempdir'
        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)

        self.hmmcopy_ref = os.path.join(ref, 'read_counter_hmm.wig')
        self.bam = os.path.join(ref, 'test.bam')
        self.output = os.path.join(self.tempdir, 'readcounter_out.txt')
        self.window_size = 1000
        self.chromosomes = ['Y']

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def run_readcounter(self):
        ReadCounter(
            self.bam,
            self.output,
            self.window_size,
            self.chromosomes,
            20,
            True).main()

    def compare_wigs(self):
        pydata = {}

        mismatches = 0
        for line in open(self.output):
            if line.startswith('data'):
                continue

            line = line.strip().split()

            chrom = line[1]
            start = int(line[2])
            stop = int(line[3])
            count = int(line[4])
            pydata[(chrom, start, stop)] = count

        for line in open(self.hmmcopy_ref):
            if line.startswith('data'):
                continue

            line = line.strip().split()
            chrom = line[1]
            start = int(line[2])
            stop = int(line[3])

            count = int(line[4])
            pycount = pydata[(chrom, start, stop)]

            if count != pycount:
                mismatches += 1

        # some reads are inexplicable missed by hmmcopy, in this case one read
        # is missed.
        self.assertLess(mismatches, 2, 'more than 1 bins have mismatches')

    def test(self):
        self.run_readcounter()
        self.compare_wigs()


class TestGCCounter(unittest.TestCase):

    def setUp(self):

        ref = 'ref'
        self.tempdir = 'tempdir'
        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)

        self.hmmcopy_ref = os.path.join(ref, 'gc_counter_hmm.wig')
        self.reference = os.path.join(ref, 'GRCh37-lite.fa')
        self.output = os.path.join(self.tempdir, 'gccounter_out.txt')
        self.window_size = 1000
        self.chromosomes = ['Y']

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def run_gccounter(self):
        GCCounter(
            self.reference,
            self.output,
            self.chromosomes,
            self.window_size).main()

    def compare_wigs(self):
        # assumes single chr in wig
        with open(self.hmmcopy_ref) as reference:
            header = reference.readline().strip().split()

            chrom = header[1].split('=')[1]
            start = header[2].split('=')[1]
            step = header[3].split('=')[1]

            data = [float(l.strip()) for l in reference]

        with open(self.output) as infile:
            header = infile.readline().strip().split()

            assert chrom == header[1].split('=')[1]
            assert start == header[2].split('=')[1]
            assert step == header[3].split('=')[1]

            indata = [float(l.strip()) for l in infile]

            self.assertEqual(data, indata, 'results did not match')

    def test(self):
        self.run_gccounter()
        self.compare_wigs()


class TestMapCounter(unittest.TestCase):

    def setUp(self):

        ref = 'ref'
        self.tempdir = 'tempdir'
        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)

        self.hmmcopy_ref = os.path.join(ref, 'map_counter_hmm.wig')
        self.reference = os.path.join(ref, 'GRCh37-lite.fa')
        self.output = os.path.join(self.tempdir, 'mapcounter_out.txt')
        self.window_size = 1000
        self.chromosomes = ['Y']

    def tearDown(self):
        pass
#         shutil.rmtree(self.tempdir)

    def run_gccounter(self):
        MapCounter(self.tempdir, self.output, self.reference, self.chromosomes,
                   'bowtie', 4, self.window_size, self.window_size, False).main()

    def compare_wigs(self):
        # assumes single chr in wig
        with open(self.hmmcopy_ref) as reference:
            header = reference.readline().strip().split()

            chrom = header[1].split('=')[1]
            start = header[2].split('=')[1]
            step = header[3].split('=')[1]

            data = [float(l.strip()) for l in reference]

        with open(self.output) as infile:
            header = infile.readline().strip().split()

            assert chrom == header[1].split('=')[1]
            assert start == header[2].split('=')[1]
            assert step == header[3].split('=')[1]

            indata = [float(l.strip()) for l in infile]

            self.assertEqual(data, indata, 'results did not match')

    def test(self):
        self.run_gccounter()
        self.compare_wigs()

if __name__ == "__main__":

    suite = unittest.TestSuite()
    loader = unittest.TestLoader()

    gc = loader.loadTestsFromTestCase(TestGCCounter)
    read = loader.loadTestsFromTestCase(TestReadCounter)
    mapc = loader.loadTestsFromTestCase(TestMapCounter)

    suite.addTests(gc)
    suite.addTests(read)
    suite.addTests(mapc)

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
