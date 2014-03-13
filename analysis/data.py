import os
import numpy
import scipy_data_fitting

class Fig4(scipy_data_fitting.Data):
    """
    Use this to load the data from Figure 4 in PhysRevLett.105.167202.

    Should not be used directly, but only subclassed.
    """

    def __init__(self, subfig):
        super().__init__()
        self.subfig = subfig
        self.name = 'fig_4' + self.subfig
        self.genfromtxt_args['delimiter'] = "\t"
        self.genfromtxt_args['skip_header'] = 1
        self.path = os.path.join('data', 'PhysRevLett.105.167202',
            'figure_4' + subfig + '.tsv')
        if subfig == 'd': self.scale = (1, 'milli')

class Fig4Parallel(Fig4):
    """
    The parallel field data from Figure 4 in PhysRevLett.105.167202.
    """

    def __init__(self, subfig):
        super().__init__(subfig)
        self.name = self.name + '_parallel'
        self.genfromtxt_args['usecols'] = (0, 1)
        if subfig == 'c': self.path = self.path.replace('.tsv', '.1.tsv')

class Fig4Antiparallel(Fig4):
    """
    The antiparallel field data from Figure 4 in PhysRevLett.105.167202.
    """

    def __init__(self, subfig):
        super().__init__(subfig)
        self.name = self.name + '_antiparallel'
        if subfig == 'c':
            self.path = self.path.replace('.tsv', '.2.tsv')
        else:
            self.genfromtxt_args['usecols'] = (0, 2)

class Fig4Difference(scipy_data_fitting.Data):
    """
    The difference of the parallel and antiparallel field data
    from Figure 4 in PhysRevLett.105.167202.
    """

    def __init__(self, subfig):
        super().__init__()
        self.subfig = subfig
        self.name = 'fig_4' + self.subfig + '_difference'
        parallel = Fig4Parallel(self.subfig)
        antiparallel = Fig4Antiparallel(self.subfig)
        self.array = numpy.array([
            parallel.array[0],
            abs(parallel.array[1] - antiparallel.array[1])
        ])

class Fig4Normalized(scipy_data_fitting.Data):
    """
    The normalized field data from Figure 4 in PhysRevLett.105.167202.

    Should not be used directly, but only subclassed.
    """

    def __init__(self, subfig, data_class):
        super().__init__()
        self.subfig = subfig
        self.name = 'fig_4' + self.subfig + '_normalized'
        self.unnormalized = data_class(self.subfig).array

        self.array = numpy.array([
            self.unnormalized[0],
            self.unnormalized[1] / max(abs(self.unnormalized[1]))
        ])

class Fig4NormalizedParallel(Fig4Normalized):
    """
    The normalized parallel field data from Figure 4 in PhysRevLett.105.167202.
    """

    def __init__(self, subfig):
        super().__init__(subfig, Fig4Parallel)
        self.name = self.name + '_parallel'

class Fig4NormalizedAntiparallel(Fig4Normalized):
    """
    The normalized antiparallel field data from Figure 4 in PhysRevLett.105.167202.
    """

    def __init__(self, subfig):
        super().__init__(subfig, Fig4Antiparallel)
        self.name = self.name + '_antiparallel'

class Fig4NormalizedDifference(Fig4Normalized):
    """
    The difference of the normalized parallel and antiparallel field data
    from Figure 4 in PhysRevLett.105.167202.
    """

    def __init__(self, subfig):
        super().__init__(subfig, Fig4Difference)
        self.name = self.name + '_difference'
