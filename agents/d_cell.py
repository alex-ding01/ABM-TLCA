import agents.farmer as farmer_class
import agents.company as company_class

class Dcell:
    """A class that contains both static and dynamic information about the simulation domain


    """
    def __init__(self, **kwargs):
        self.nFarmAgent = 0
        self.nCompanyAgent = 0
        self.FarmerAgents = []
        self.CompanyAgents = []

    def add_agent(self, agent_struct):
        """Adds a new agent to a list
        """
        agent_type = type(agent_struct).__name__

        if agent_type == company_class.Company.__name__:
            self.nCompanyAgent += 1
            self.CompanyAgents.append(agent_struct)

        if agent_type == farmer_class.Farmer.__name__:
        # if agent_type == 'Farmer':  
            self.nFarmAgent += 1
            self.FarmerAgents.append(agent_struct)

        return agent_type
