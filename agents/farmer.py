

class Farmer:

    def __init__(self, **kwargs):

        self.FarmSize = kwargs.get('FarmSize')
        self.CropType = kwargs.get('CropType')
        self.AgentType = kwargs.get('AgentType')

        self.SizeAlpha = kwargs.get('SizeAlpha')
        self.SizeBeta = kwargs.get('SizeBeta')


    def UpdateFamiliarMisc2chp(self, NewFamiliarMiscAlpha2chp, NewFamiliarMiscBeta2chp):
        """ Updates the farmers likelihood of switching their crop such that they are less likely to switch as they age
         """
        self.FamiliarMiscAlpha2chp = FamiliarMiscAlpha2chp
        self.FamiliarMiscBeta2chp = FamiliarMiscBeta2chp

    def UpdateFamiliarMisc2bm(self, NewFamiliarMiscAlpha2bm, NewFamiliarMiscBeta2bm):
        """ Updates the farmers likelihood of switching their crop such that they are less likely to switch as they age
         """
        self.FamiliarMiscAlpha2bm = NewFamiliarMiscAlpha2bm
        self.FamiliarMiscBeta2bm = NewFamiliarMiscBeta2bm

