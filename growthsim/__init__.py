# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 14:39:12 2016

@author: carolc24
"""

import simplesbml
import libsbml
import roadrunner

class strain(simplesbml.sbmlModel):
    def __init__(self):
        reader = libsbml.SBMLReader();
        doc = reader.readSBMLFromFile("cellmodel.xml");
        self.document = doc;
        self.model = doc.getModel();
        
    def addProtein(self, protein_id, n=300, txn_rate="we*a/(thetax + a)", modifiers=["a"], deg_rate="0"):
        self.addSpecies('[' + protein_id + ']', 0, "cell");
        self.addSpecies('[m' + protein_id + ']', 0, "cell");
        self.addSpecies('[rm' + protein_id + ']', 0, "cell");
        self.getModel().getAssignmentRule("ttrate").getMath().getChild(0).prependChild( \
                libsbml.parseL3Formula("rm" + protein_id));
        self.getModel().getAssignmentRule("fr").getMath().getChild(0).getChild(1). \
                prependChild(libsbml.parseL3Formula("rm" + protein_id));
        self.getModel().getAssignmentRule("fr").getMath().getChild(1).getChild(0). \
                getChild(1).prependChild(libsbml.parseL3Formula("rm" + protein_id));
        self.getModel().getAssignmentRule("fr").getMath().getChild(1).getChild(1). \
                getChild(1).prependChild(libsbml.parseL3Formula(protein_id));
    
        self.addParameter("n" + protein_id, n, "cell");
        v = self.addReaction([], ['m' + protein_id], txn_rate);
        for s in modifiers:
            m = v.createModifier();
            m.setSpecies(s);
        self.addReaction(['m' + protein_id, "r"], ['rm' + protein_id], "kb * r * m" + protein_id);
        self.addReaction(['rm' + protein_id], ['r', 'm' + protein_id], "ku * rm" + protein_id);
        self.addReaction(['rm' + protein_id], ['r', 'm' + protein_id, protein_id], \
                    "gamma * rm" + protein_id + " / n" + protein_id);
        self.addReaction(['m' + protein_id], [], 'lam * m' + protein_id);
        self.addReaction(['rm' + protein_id], [], 'lam * rm' + protein_id);
        self.addReaction([protein_id], [], 'lam * ' + protein_id);
        self.addReaction(['m' + protein_id], [], 'dm * m' + protein_id);
        self.addReaction([protein_id], [], deg_rate + " * " + protein_id);

    def change_ids(self, elem_list, tag):
        for elem in elem_list:
            old_id = elem.getId();
            elem.setId(old_id + "_" + tag);        
          
    def change_spec(self, elem_list, tag):
        for elem in elem_list:
            old_id = elem.getSpecies();
            elem.setSpecies(old_id + "_" + tag);
            
    def change_rule(self, elem_list, tag):
        for elem in elem_list:
            old_id = elem.getVariable();
            elem.setVariable(old_id + "_" + tag);
            math = elem.getMath();
            self.parse_rule(math, tag);
            
    def change_reac(self, elem_list, tag):
        for elem in elem_list:
            old_id = elem.getId();
            elem.setId(old_id + "_" + tag);
            math = elem.getKineticLaw().getMath();
            self.parse_rule(math, tag);
            self.change_spec(elem.getListOfReactants(), tag);
            self.change_spec(elem.getListOfProducts(), tag);
            self.change_spec(elem.getListOfModifiers(), tag);
            
    def parse_rule(self, math, tag):
        if math:
            name = math.getName();
            if name:
                math.setName(name + "_" + tag);
            for c in range(math.getNumChildren()):
                self.parse_rule(math.getChild(c), tag);
                
    def process_model(self, tag):
        rr_init = roadrunner.RoadRunner(self.toSBML());
        result = rr_init.simulate(0,5000,101);
        for i in range(len(self.getListOfSpecies())):
            s = self.getSpecies(i);
            s.setInitialConcentration(result[100,i]);
        self.addSpecies("N", 1, "cell");
        self.getSpecies("et").setInitialConcentration(1.0);
        self.getParameter("Kp").setValue(1260./(3*10**8));
        self.getParameter("ns").setValue(100);
        self.change_ids(self.getListOfSpecies(), tag);
        self.change_ids(self.getListOfParameters(), tag);
        self.change_reac(self.getListOfReactions(), tag);
        self.change_rule(self.getListOfRules(), tag);
    
    def clone(self):
        mod = self.model.clone();
        new_strain = strain();
        new_strain.document.setModel(mod);
        new_strain.model = mod;     
        return new_strain;

class population(simplesbml.sbmlModel):    
    def __init__(self, tag_dict):
        self.strains = dict();
        mods = tag_dict.values();
        tags = tag_dict.keys();
        
        simplesbml.sbmlModel.__init__(self);
        self.addCompartment(1.0, "cell");
        
        self.addSpecies('s', 1e8, "cell");
        self.addParameter("dn", 0.03);
        self.addParameter("k_in", 0.03*1e8);
        self.addReaction([], ["s"], "k_in");
        self.addReaction(["s"], [], "dn*s");
        
        for i in range(len(mods)):
            self.addStrain(mods[i], tags[i]);
        
    def addStrain(self, modn, tag):
        if (self.strains.has_key(tag)):
            raise SystemExit("Error: this tag is already used");
            
        self.strains.update({tag:modn});
        modn.process_model(tag);

        self.getListOfSpecies().appendFrom(modn.getListOfSpecies());
        self.getListOfParameters().appendFrom(modn.getListOfParameters());
        self.getListOfReactions().appendFrom(modn.getListOfReactions());
        self.getListOfRules().appendFrom(modn.getListOfRules());

        kl = "N_{0}*et_{0}*vt_{0}*s/(Kt_{0}+s)".format(tag);
        v = self.addReaction(["s"], [], kl);
        m = v.createModifier();
        m.setSpecies("N_" + tag);
        m = v.createModifier();
        m.setSpecies("et_" + tag);
        gro = "lam_{0}*N_{0} - dn*N_{0}".format(tag);
        self.addRateRule("N_" + tag, gro);
    
        math = self.getReaction("r_0_" + tag).getKineticLaw().getMath();
        math.getLeftChild().getRightChild().setName("s");
        math.getRightChild().getRightChild().setName("s");
        m = self.getReaction("r_0_" + tag).createModifier();
        m.setSpecies("s");
        
    def simulate(self, tstart, tend, nsteps):
        rr = roadrunner.RoadRunner(self.toSBML());
        pop_tags = sorted(self.strains.keys());
        selections = ['time'];
        for i in range(len(self.strains)):
            selections.append('N_' + pop_tags[i]);
        rr.timeCourseSelections = selections;
        result = rr.simulate(tstart, tend, nsteps);
        return result;