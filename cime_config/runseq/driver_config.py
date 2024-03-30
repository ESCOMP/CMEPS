#!/usr/bin/env python3

# Inherit from the dictionary class
class DriverConfig(dict):

    ###############################################
    def __init__(self, case, coupling_times):
    ###############################################
        # this initializes the dictionary
        super(DriverConfig,self).__init__()

        self['atm'] = self.__compute_atm(case, coupling_times)
        self['glc'] = self.__compute_glc(case, coupling_times)
        self['ice'] = self.__compute_ice(case, coupling_times)
        self['lnd'] = self.__compute_lnd(case, coupling_times)
        self['ocn'] = self.__compute_ocn(case, coupling_times)
        self['rof'] = self.__compute_rof(case, coupling_times)
        self['wav'] = self.__compute_wav(case, coupling_times)

    ###############################################
    def __compute_atm(self, case, coupling_times):
    ###############################################

        comp_atm  = case.get_value("COMP_ATM")

        run_atm = True
        med_to_atm = True
        if (comp_atm == 'satm'):
            run_atm = False
            med_to_atm = False
        elif (comp_atm == 'datm'):
            # TODO: check of data model prognostic flag is on - this is a new xml variable
            # If the prognostic flag is on, then should set med_to_atm to True
            if_prognostic = False
            med_to_atm = if_prognostic

        # TODO: need to put in special logical for adiabatic mode - where the atmosphere does
        # does not really send anything to the mediator
        # This will be the case if (case.get_value('COMP_OCN') == 'socn'):
        # In this case - a special run sequence should be written

        return (run_atm, med_to_atm, coupling_times["atm_cpl_dt"])

    ###############################################
    def __compute_glc(self, case, coupling_times):
    ###############################################

        # In the mediator the glc_avg_period will be set as an alarm
        # on the on the prep_glc_clock. When this alarm rings - the
        # averaging will be done.

        comp_glc = case.get_value("COMP_GLC")

        run_glc = True
        med_to_glc = True
        if (comp_glc == 'sglc'):
            run_glc = False
            med_to_glc = False
        elif (comp_glc == 'cism'):
            if not case.get_value("CISM_EVOLVE"):
                run_glc = False

        # If CISM is not evolving only get data back from cism at the initial time
        # However will still need to call the exchange at the end if the stop_option
        # is nsteps or days - or otherwise just every ndays
        # Note that nsteps is the minimum component coupling time
        if (comp_glc == 'cism'):
            glc_coupling_time = coupling_times["glc_cpl_dt"]
            if not case.get_value("CISM_EVOLVE"):
                stop_option = case.get_value('STOP_OPTION')
                stop_n = case.get_value('STOP_N')
                if stop_option == 'nyears':
                    glc_coupling_time = coupling_times["glc_cpl_dt"]
                elif stop_option == 'nsteps':
                    glc_coupling_time = stop_n * coupling_times["glc_cpl_dt"]
                elif stop_option == 'ndays':
                    glc_coupling_time = stop_n * 86400
                else:
                    glc_coupling_time = 86400
        elif (comp_glc == 'dglc'):
            glc_coupling_time = coupling_times["glc_cpl_dt"]
            stop_option = case.get_value('STOP_OPTION')
            stop_n = case.get_value('STOP_N')
            if stop_option == 'nsteps':
                glc_coupling_time = stop_n*coupling_times["atm_cpl_dt"]
        elif (comp_glc == 'xglc'):
            glc_coupling_time = coupling_times["glc_cpl_dt"]
        else:
            glc_coupling_time = 0

        return (run_glc, med_to_glc, glc_coupling_time)

    ###############################################
    def __compute_ice(self, case, coupling_times):
    ###############################################

        comp_ice  = case.get_value("COMP_ICE")

        run_ice = True
        med_to_ice = True # this is the case for dice
        if (comp_ice == 'sice'):
            run_ice = False
            med_to_ice = False

        return (run_ice, med_to_ice, coupling_times["ice_cpl_dt"])

    ###############################################
    def __compute_lnd(self, case, coupling_times):
    ###############################################

        comp_lnd  = case.get_value("COMP_LND")

        run_lnd = True
        med_to_lnd = True
        if (comp_lnd == 'slnd'):
            run_lnd = False
            med_to_lnd = False
        elif (comp_lnd == 'dlnd'):
            # TODO: check of data model prognostic flag is on - this is a new xml variable
            # If the prognostic flag is on, then should set med_to_lnd to True
            if_prognostic = False
            med_to_lnd = if_prognostic

        return (run_lnd, med_to_lnd, coupling_times["lnd_cpl_dt"])

    ###############################################
    def __compute_ocn(self, case, coupling_times):
    ###############################################

        comp_ocn = case.get_value("COMP_OCN")

        run_ocn = True
        med_to_ocn = True
        if (comp_ocn == 'socn'):
            run_ocn = False
            med_to_ocn = False
        elif (comp_ocn == 'docn'):
            # TODO: check of data model prognostic flag is on - this is a new xml variable
            # If the prognostic flag is on, then should set med_to_wav to True
            docn_mode = case.get_value("DOCN_MODE")
            docn_import_fields = case.get_value("DOCN_IMPORT_FIELDS")
            med_to_ocn = ('som' in docn_mode or 'interannual' in docn_mode or docn_import_fields != 'none')

        return (run_ocn, med_to_ocn, coupling_times["ocn_cpl_dt"])

    ###############################################
    def __compute_rof(self, case, coupling_times):
    ###############################################

        comp_rof  = case.get_value("COMP_ROF")

        run_rof = True
        med_to_rof = True
        if (comp_rof == 'srof'):
            run_rof = False
            med_to_rof = False
        elif (comp_rof == 'drof'):
            # TODO: check of data model prognostic flag is on - this is a new xml variable
            # If the prognostic flag is on, then should set med_to_rof to True
            if_prognostic = False
            med_to_rof = if_prognostic
        else:
            # this is active runoff - determine if the mode or the grid is null - and in that case
            # remove all interactions with rof from the run sequence
            if ((case.get_value("COMP_ROF") == 'mosart' and case.get_value("MOSART_MODE") == 'NULL') or
                (case.get_value("COMP_ROF") == 'rtm' and case.get_value("RTM_MODE") == 'NULL') or
                (case.get_value("ROF_GRID") == 'null')):
                run_rof = False
                med_to_rof = False

        return (run_rof, med_to_rof, coupling_times["rof_cpl_dt"])

    ###############################################
    def __compute_wav(self, case, coupling_times):
    ###############################################

        comp_wav  = case.get_value("COMP_WAV")

        run_wav = True
        med_to_wav = True
        if (comp_wav == 'swav'):
            run_wav = False
            med_to_wav = False
        elif (comp_wav == 'dwav'):
            # TODO: check of data model prognostic flag is on - this is a new xml variable
            # If the prognostic flag is on, then should set med_to_wav to True
            if_prognostic = False
            med_to_wav = if_prognostic

        return (run_wav, med_to_wav, coupling_times["wav_cpl_dt"])
