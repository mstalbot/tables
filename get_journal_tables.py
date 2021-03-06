from astropy.table import Table
from bs4 import BeautifulSoup
from os.path import join, exists
from requests import Session
from astropy.io import ascii
from collections import OrderedDict
from astropy.coordinates import SkyCoord
import numpy as np
import json
import pandas as pd
import re
from time import sleep
import csv
from astropy.io import fits
import xml.dom.minidom
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import hms_tuple, dms_tuple
from astropy import units

class Journal_tables():
 
    def __init__(self, get_bibcodes = False, get_online_tables = False, user_name = 'guest', base_directory = 'PATH_TO_lens_surveys', headers = {'user-agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_6) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/13.1.2 Safari/605.1.15'}, user = {'MNRAS':'','IOP':'','A&A':''}, password = {'MNRAS':'','IOP':'','A&A':''}, mld_auth = {'user':'','password':''}, start = 0, end = 99999, redo_pandas = False, rescan_online = False, slow_down_seconds_after_requests = 5, inspect = False, redo_inspection=False, load_processed_data = True, prepare_to_post_lenses_to_MLD2 = False):
    
        """This program retrieves tables from journals typically used in Astronomy. Data is saved as a pandas or ascii version of a JSON file.
        
        Usage:
        get_bibcodes: Reads in bibcodes from a json file
        get_online_tables: Download tables from journals
        headers: Headers from your site. The provided may be sufficient.
        user: User logins for journals
        password: Passwords for journals
        base_directory: Where json files of journal tables are saved and loaded from.
        start: Where to start process in bibcode list
        end: Where to end process in bibcode list
        rescan_online: Rescan online data (usefull in case requests were blocked for overuse.
        redo_pandas: Redo pandas tables for touch ups if a change in that process is made
        slow_down_seconds_after_requests: Crawl speed for internet requests
        inspect: Inspect tables not yet inspected
        redo_inspection: Redo inspection
        mld_auth: Masterlens database login in format of {'user':'Your user login', 'password':'Your password'}.
        load_processed_data: Load previous processed data and relavent MLD database information
        """ 
        
        #Set class globals
        self.headers = headers
        self.get_bibcodes = get_bibcodes
        self.base_directory = base_directory
        self.headers = headers
        self.get_online_tables = get_online_tables
        self.base_links = {'ADS': 'https://ui.adsabs.harvard.edu', 'IOP': 'https://iopscience.iop.org', 'MNRAS': 'https://academic.oup.com', 'A&A':'https://www.aanda.org'}
        self.journal_id_converter = {"0":"", "44":"Astrofizicheskie Issledovaniia Izvestiya Spetsial'noj Astrofizicheskoj Observatorii", "37":"arXiv e-prints", "11":"Astronomy and Astrophysics", "10":"The Astronomical Journal", "1":"The Astrophysical Journal", "8":"Astrophys. J. Lett.", "42":"The Astrophysical Journal Supplement Series", "45":"Bulletin d'information du telescope Canada-France-Hawaii", "2":"Monthly Notices of the Royal Astronomical Society", "39":"Nat", "38":"Nuc. Phy. B Proc. S.", "40":"Publications of the Astronomical Society of Japan", "43":"Revista Mexicana de Astronomia y Astrofisica", "41":"Science"}
        self.journal_id_converter_inverted = {self.journal_id_converter[key]:key for key in self.journal_id_converter}
        self.start = start
        self.end = end
        self.redo_pandas = redo_pandas
        self.rescan_online = rescan_online
        self.inspect = inspect
        self.redo_inspection = redo_inspection
        self.load_processed_data = load_processed_data
        self.prepare_to_post_lenses_to_MLD2 = prepare_to_post_lenses_to_MLD2
        self.slow_down_seconds_after_requests = slow_down_seconds_after_requests
        self.fix_messed_up_tables = {'2016ApJ...826..112S': '\n\n\n\n\n\n\n\n\t'}
        self.table_action = {'^':'Word to recognize name is of lens and NOT source', '~':'Name,Ra,Dec of cluster or group lens','+':'Cluster Sources Table', 's': 'Skip cause table not important', 'y': 'Table correctly loaded', 'h': 'Append first table row to header', 'f': 'Trim footer', 'p': 'Warn problem', 'm': 'Message on table quality', 'c': 'Check system name since may be unique', '?':'Table might not be necessary', 't':'Detection table', 'i':'Follow up table', '-': 'Double check diction', 'd':'Skipping since duplicate', 'z':'Problem but check. Currently skipping', 'r':'Repeats format horizontally', 'n':'Discovery', 'k':'Lens Kind', 'g': 'Candidate table represents a grade', 'q':'Table columns are horizontal', '=':'Word to recognize candidate in column'}
        self.table_inspection_id = {'sn': 'Source names', 's': 'System Name', 'dd': 'Discovery Date', 'an': 'Alternate Name(s)', 'ndp': 'Number of Discovery Programs', 'rh': 'RA in Hours', 'rf': 'RA in Hours:Min:Sec', 'rpw': 'RA in Degrees:Min:Sec', 'rhp': 'RA (Hours part)', 'rmp': 'RA (Mins part)', 'rsp': 'RA (Secs part)', 'r': 'RA [°]', 'dp': 'Dec (+/-) Degree:Min:Sec', 'ddp': 'Dec (Degree part)', 'dmp': 'Dec (Arcmin part)', 'dsp': 'Dec (Arcsec part)', 'd': 'Dec [°]', 'lg': 'Lens Grade', 'er': 'Einstein_R ["]', 'erq': 'Einstein_R ["] quality', 'ere':'Einstein_R ["] error', 'zl': 'z_Lens', 'zlq': 'z_Lens quality', 'zle':'z_Lens error', 'zs': 'z_Source(s)', 'zsq': 'z_Source quality', 'zse':'z_Source error', 'svd': 'Stellar velocity disp', 'des': 'Description', 'ref': 'References', 'el': 'External Links', 'dg': 'Detection score or grade', 'p':'Plate', 'mjd':'MJD', 'f':'Fiberid', 'pos':'Position', 'r-d':'RA-Dec (Degrees)', 'p2':'Plate2', 'mjd2': 'MJD2', 'f2':'Fiberid2', 'p3':'Plate3', 'mjd3': 'MJD3', 'f3': 'Fiberid3', 'zl2': 'z_Lens2', 'zl3': 'z_Lens3', 'es':'Has external link for SDSS','eads':'Has external link for ADS','en':'Has external link for NED','eapod':'Has external link for APOD','svde':'Stellar velocity disp error','dc':'Discovery count', 'kr':'killreferences', 'ar':'addreferences', 'rk':'referencestokill[]', 'ra':'referencestoadd[]', 'lt': 'Lens type', 'ic':'Image count', 'pmf':'Plate-MJD-Fiberid', 'c':'Condition'}
        
        #User MLD authorization data, urls, and plain english to HTML entry conversion metadata
        self.mld_auth = mld_auth
        
        self.masterlens_url = 'http://test.masterlens.org'
        self.masterlens_login_url = join(self.masterlens_url, 'member.php?action=login&')
        self.masterlens_form_extension = 'search.php?'
        self.masterlens_form_citation = 'citation.php?'
        self.masterlens_form_mode = 'mode=add&'
        self.masterlens_form_add_a_lens = join(self.masterlens_url, self.masterlens_form_extension, self.masterlens_form_mode)
        self.masterlens_form_add_a_paper = join(self.masterlens_url, self.masterlens_form_citation, self.masterlens_form_mode)
        self.masterlens_phrases_to_input_converter = {"Submit":"inputaction", "Has external link for SDSS":"query_has_sdss", "Has external link for ADS":"query_has_adsabs", "Has external link for NED":"query_has_ned", "Has external link for APOD":"query_has_apod", "APOD link":"query_apod_link", "Lens Grade":"query_lensgrade", 'Einstein_R ["]':"query_theta_e", 'Einstein_R ["] error':"query_theta_e_err", 'Einstein_R ["] quality':"query_theta_e_quality", "z_Lens":"query_z_lens", "z_Lens error":"query_z_lens_err", "z_Lens quality":"query_z_lens_quality", "z_Source(s)":"query_z_source", "z_Source error":"query_z_source_err", "z_Source quality":"query_z_source_quality", "Stellar velocity disp":"query_vdisp", "Stellar velocity disp error":"query_vdisp_err", "Discovery count":"query_discovery_count", "query_discoveryID_num1":"0", "Lens type":"query_kindID", "Description":"query_description", "Standard Name":"query_system_name", "Discovery Date":"query_discovery_date", "System Name":"query_alternate_name", "RA (Hours part)":"query_ra_hrs", "RA (Mins part)":"query_ra_mins", "RA (Secs part)":"query_ra_secs", "RA [°]":"query_ra_coord", "Dec (Degree part)":"query_dec_degrees", "Dec (Arcmin part)":"query_dec_arcmin", "Dec (Arcsec part)":"query_dec_arcsec", "Dec [°]":"query_dec_coord", "killref":"killreferences", "addreferences": "addreferences", "killreferences":"referencestokill[]", "reference":"referencestoadd[]", "Kind_MLD_ID": "query_kindID", "Discovery_MLD_ID": "query_discoveryID_num1"}
        
        #This enables diction keys for the inspection process, including MLD post keys that are not loaded in an MLD xml file.
        self.quick_column_letter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','x','y','z','aa','bb','cc','dd','ee','ff','gg','hh','ii','jj','kk','ll','mm','nn','oo','pp','qq','rr','ss','tt','uu','vv','xx','yy','zz']
        
        self.z_type_id = {"":"0", "spectroscopic":"1", "photometric":"2"}
        self.er_quality_id = {"":"0", "SIE model":"1", "1/2 image separation":"2", "Reference redshift":"3"}        
        self.z_type_id_inverted = {self.z_type_id[key]:key for key in self.z_type_id}
        self.er_quality_id_inverted = {self.er_quality_id[key]:key for key in self.er_quality_id}

        #Set global overview and tables dictionary to save and load data from
        self.ads_scrapped_data = {}
        self.ads_scrapped_tables = {}
        
        #Automatically run if True
        if self.get_bibcodes: self.load_query_bibcodes()
        if self.load_processed_data: self.load_saved_data()
        if self.get_online_tables:
            self.set_sessions('open', user, password)
            self.set_papers_overview()
            self.set_tables()
            self.set_sessions('close')
            self.report_overal_stats()
        if self.inspect: self.inspect_tables()
        if self.prepare_to_post_lenses_to_MLD2: self.set_system_data()
    
    
    def load_query_bibcodes(self):
        """Load bibliography codes from a json file on disk"""
        
        with open(join(self.base_directory, 'resources', 'bibcodes_leonidas.txt'), 'r') as json_file:
            self.bibcodes = json.load(json_file)['bibcodes']
            
                
    def set_sessions(self, open_or_close='open', user={'MNRAS':'','IOP':'','A&A':''}, password={'MNRAS':'','IOP':'','A&A':''}):
        """Set sessions to web pages with attempt to login"""
        
        if open_or_close == 'open':
            self.sessions = {'ADS': Session()}
            for base in self.base_links.keys():
                if 'ADS' not in base:
                    self.sessions[base] = Session()
                    self.sessions[base].auth = (user[base], password[base])
                    self.sessions[base].post(self.base_links[base], headers = self.headers)
        elif open_or_close == 'close':
            for base in self.base_links.keys():
                self.sessions[base].close()
                del self.sessions[base]
                
                
    def load_saved_data(self):
        """Load saved data"""

        #Load in reference data
        self.load_MLD_surveys_ids()
        self.load_MLD_kinds_ids()
        self.load_MLD_references()
        self.load_MLD_lenses()
        self.load_sugohi()
        self.load_silo_eboss()
        self.load_links()
        
        #Load in processed data
        for self.query in self.bibcodes:
            self.load_saved_overview()
            self.load_saved_tables()
            
    def get_inspection_status(self):
        """Report completed table counts and skipped tables"""
        
        table_count=0
        skipped = 0
        
        #Loop through papers
        for index, self.query in enumerate(self.bibcodes[self.start:self.end]):
            print('INDEX>>>', index)
            
            #Load saved table if not already done.
            self.load_saved_overview()
            self.load_saved_tables()
            
            #Loop through tables to find inspection status
            if ('Status' in self.ads_scrapped_data[self.query] and self.ads_scrapped_data[self.query]['Status'] == 'Complete'):
                for key in sorted(self.ads_scrapped_data[self.query]['Table meta data'].keys()):
                    if 'Pandas format' in self.ads_scrapped_data[self.query]['Table meta data'][key]:
                        for key2 in sorted(self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'].keys()):
                            if 'Inspection' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]:
                                for ikey in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']:
                                    if 'Skip' in ikey or 'skipping' in ikey or 'Skipping' in ikey:
                                        skipped += 1
                                        break
                            table_count+=1
            print('Total tables inspected:', table_count, ', Tables considered relavent:', table_count-skipped)
                            
                
    def inspect_tables(self):
        """Inspect tables"""
        
        to_scan = 0
        #Loop through papers
        for index, self.query in enumerate(self.bibcodes[self.start:self.end]):
            print('INDEX>>>', index)
            
            #Load saved table if not already done.
            self.load_saved_overview()
            self.load_saved_tables()
            
            #Check if inspection of the tables is needed
            if ('Status' in self.ads_scrapped_data[self.query] and self.ads_scrapped_data[self.query]['Status'] == 'Complete') and ('Inspection status' not in self.ads_scrapped_data[self.query] or self.ads_scrapped_data[self.query]['Inspection status'] != 'Complete' or self.redo_inspection):
                
                #Display table and overview data.
                if 'Paper Overview' in self.ads_scrapped_data[self.query]:
                    try: print('\n=========================================Paper Title=========================================\n', self.ads_scrapped_data[self.query]['Paper Overview']['Title'], '\n=========================================Paper Authors=========================================\n', self.ads_scrapped_data[self.query]['Paper Overview']['authors'], '\n=========================================Paper Abstract=========================================\n', self.ads_scrapped_data[self.query]['Paper Overview']['abstract'], '\n', 'Bibcode:', self.query)
                    except: print('\n=========================================Bibcode: %s=========================================\n'%self.query)
                else: print('\n=========================================Bibcode: %s=========================================\n'%self.query)
                
                unbroken = True 
                for key in sorted(self.ads_scrapped_data[self.query]['Table meta data'].keys()):
                    if self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] == 'Complete':
                        print('Paper link:', self.ads_scrapped_data[self.query]['Publisher link via ADS gateway'])
                        print('Table link:', self.ads_scrapped_data[self.query]['Table meta data'][key]['Link'])
                        for key2 in sorted(self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'].keys()):
                            #Display table information
                            redo = True
                            while redo:
                                if self.query not in self.ads_scrapped_data[self.query]['Table meta data'][key]['Link']:
                                    print('\n--------------------------------------------------Extra info--------------------------------------------------\n')
                                    for ekey in self.ads_scrapped_data[self.query]['Table meta data'][key]['Table captions or footers'].keys(): print('%s:'%ekey, self.ads_scrapped_data[self.query]['Table meta data'][key]['Table captions or footers'][ekey])
                                print('\n----------------------------------------------------Table-----------------------------------------------------\n', self.ads_scrapped_tables[self.query][key][key2])
                                #Will skip table if inspection already done and re-inspect not stated.
                                if 'Inspection' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2] and not self.redo_inspection:
                                    print('\n\n>>>Already inspected. Skipping\n\n')
                                    redo = False
                                    continue
                                    
                                if 'Inspection' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2] and ('Skip cause table not important' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'] or 'Skipping since duplicate' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'] or 'Problem but check. Currently skipping' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']):
                                    print('Skipping cause noted to skip')
                                    redo = False
                                    continue
                                    
                                '''if 'Inspection' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2] and ('Skip cause table not important' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'] or 'Problem but check. Currently skipping' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']):
                                    pass
                                else:    
                                    print('Skipping cause noted to skip AS COMPLETE')
                                    redo = False
                                    continue'''
                                
                                '''if 'Inspection' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2] and 'Message on table quality' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'] and 'source' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']['Message on table quality']:
                                    pass
                                else:
                                    print('Skipping')
                                    redo = False
                                    continue'''
                                     
                                to_scan+=1
                                #Process inspection keys for inspection display and recording.
                                columns = OrderedDict()
                                for index, mkey in enumerate(self.ads_scrapped_tables[self.query][key][key2].keys()): columns[self.quick_column_letter[index]] = mkey
                                print('>>>TABLE COLUMN QUICK KEY\n%s\n'%columns)
                                
                                
                                try: print('Last inspection save>>>', self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection'])
                                except: pass
                                #Record table quality and relavence
                                notes = {}
                                if input('Continue?') == '':
                                    redo = False
                                    continue
                                for action in input('TYPE OR ONE FROM THE FOLLOWING: %s'%self.table_action):
                                    if action not in ['f','p','m','-', 'r', 'n', 'k', 'g', '=', '^', '~']: notes[self.table_action[action]] = ''
                                    else: 
                                        comment = input(('\nDiscoverer Name from list: %s'%self.discovery_id_inverted) if 'n' in action else ('\nLens type from: %s'%self.lens_type_id_inverted) if 'k' in action else ('Conditional statement in "Condition" column to record row (add "!" to front if to reject)') if '=' in action else ('Conditional statement in "Source name" column to recognize if lens name') if '^' in action else 'Name,RA,DEC of cluster lens' if '~' in action else '%s: '%self.table_action[action])
                                        if comment == '': comment = input('Hmmm...if not specified, please state a custom identifier')
                                        notes[self.table_action[action]] = comment
                                self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection'] = {'Notes': notes}
                                
                                if 'Table columns are horizontal' in notes:
                                    columns = OrderedDict()
                                    columns['>'] = 'pdname' 
                                    for index, mkey in enumerate(self.ads_scrapped_tables[self.query][key][key2].T.keys()): columns[self.quick_column_letter[index]] = mkey
                                    print('Table rotated\n', self.ads_scrapped_tables[self.query][key][key2].T)
                                #Record column info if table not specified to be skipped.
                                if self.table_action['s'] not in notes and self.table_action['d'] not in notes and self.table_action['z'] not in notes:
                                
                                    #Initialize a column name map
                                    map = {}
                                    
                                    #Print Table key to quick inspection key
                                    print('----------------------------------------------------------------------------------------------------------\n')
                                    print('>>>TABLE COLUMN QUICK KEY\n%s\n'%columns)
                                    print('>>>Potential column names:\n', ', '.join(['%s: %s'%(mkey, self.table_inspection_id[str(mkey)]) for mkey in self.table_inspection_id.keys()]))
                                    print('----------------------------------------------------------------------------------------------------------\n')
                                    
                                    #Record user quick key input
                                    for column_id in input('Identify columns by "number:column total key" separated by ",":\n').split(','): 
                                        map[self.table_inspection_id[column_id.split(':')[0]]] = str(columns[column_id.split(':')[1]])
                                        if 'z_Lens' in self.table_inspection_id[column_id.split(':')[0]] or 'z_Source(s)' in self.table_inspection_id[column_id.split(':')[0]]:
                                            method = input('Choose zlens method type: %s'%self.z_type_id_inverted)
                                            if method != '': method = self.z_type_id_inverted[method]
                                            error = input('Choose zlens error column: %s'%columns)
                                            if error != '': error = columns[error]
                                            print(method, error, column_id.split(':'))
                                            print(column_id.split(':')[0])
                                            print(self.table_inspection_id[column_id.split(':')[0]])
                                            map[self.table_inspection_id[column_id.split(':')[0]]] += ('@' + str(method) + '@' + str(error))
                                        if 'Einstein_R ["]' in self.table_inspection_id[column_id.split(':')[0]]:
                                            method = input('Choose method type: %s'%self.er_quality_id_inverted)
                                            if method != '': method = self.er_quality_id_inverted[method]
                                            else:
                                                method = input('Is ER spec a total separation (i.e. not 1/2 sep)?')
                                                if method == 'y': method = 'Images separation'
                                                else: method = ''
                                            error = input('Choose ER error column: %s'%columns)
                                            if error != '': error = columns[error]
                                            map[self.table_inspection_id[column_id.split(':')[0]]] += ('@' + str(method) + '@' + str(error))
                                        if 'Stellar velocity disp' in self.table_inspection_id[column_id.split(':')[0]]:
                                            error = input('Choose vdisp error column: %s'%columns)
                                            if error != '': error = columns[error]
                                            map[self.table_inspection_id[column_id.split(':')[0]]] += ('@' + '' + '@' + str(error))
                                        if 'Lens type' in self.table_inspection_id[column_id.split(':')[0]]:
                                            method = input('Choose lens type method type: %s'%self.lens_type_id_inverted)
                                            if method != '': method = self.lens_type_id_inverted[method]
                                            map[self.table_inspection_id[column_id.split(':')[0]]] += ('@' + str(method) + '@' + '')
                                        #Records a discoverer name. NOTE table may not be standardized to MLD post form!!!!
                                        if self.table_inspection_id[column_id.split(':')[0]] == 'Discovery':
                                            method = input('Choose discovery method type: %s'%self.table_inspection_id_inverted)
                                            if method != '': method = self.table_inspection_id_inverted[method]
                                            else: method = input('Discovery not in MLD. Please set a acronym for now.')
                                            map[self.table_inspection_id[column_id.split(':')[0]]] += ('@' + str(method) + '@' + '')
                                        
                                    self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Table map to MasterLens database'] = map
                                print('\n>>>>>>HERE if your output:', self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection'],'\n')
                                #Save per table inspection
                                confirm = input('Please confirm correct (y=yes, n=redo)')
                                if confirm == 'y':
                                    self.save_overview()
                                    redo = False
                                else: redo = True
                                
                    else: unbroken = False
                    
                #Save if inspection of tables in paper is complete
                if unbroken: self.ads_scrapped_data[self.query]['Inspection status'] = 'Complete'
                else: self.ads_scrapped_data[self.query]['Inspection status'] = 'Incomplete'
                self.save_overview()
        print('Scanned', to_scan)

    def set_tables(self):
        """Download tables from journals"""
        
        for index, self.query in enumerate(self.bibcodes[self.start:self.end]):
            print('---------------------Run %s Running Bibcode: %s-------------------'%(index + self.start, self.query))
            
            #Load previously saved data
            self.load_saved_overview()
            self.load_saved_tables()
            
            #Determine if paper is complete and should be skipped
            if 'Status' in self.ads_scrapped_data[self.query] and self.ads_scrapped_data[self.query]['Status'] in ['Complete', 'Skipped'] and not (self.rescan_online or self.redo_pandas):
                print('>Skipping Bibcode: %s since return command is: %s'%(self.query, self.ads_scrapped_data[self.query]['Status']))
                continue
            elif 'Status' not in self.ads_scrapped_data[self.query]: self.ads_scrapped_data[self.query]['Status'] = 'Start'
            
            #Used to note if a change has been made and should be saved
            self.changed = False
            
            #Copy paper overview info from RIS output file, which file was generated from Leonidas paper list
            if self.rescan_online or 'Paper Overview' not in self.ads_scrapped_data[self.query]:
                self.ads_scrapped_data[self.query]['Paper Overview'] = self.overview_list[self.query]
                self.changed = True
                           
            #Set ADS link to journal based on ADS standard gateway format
            gateway = join(self.base_links['ADS'], 'link_gateway', self.query, 'PUB_HTML' if 'arXiv' not in self.query else 'EPRINT_HTML')
            
            #Set journal acronym from bibcode formated to have the year stated first and pages, ect separated by "."
            self.journal = self.query.split('.')[0][4:]
            
            #Set key for base requests link used
            self.base = 'IOP' if self.journal in ['AJ', 'ApJ', 'ApJS', 'NAAS'] else self.journal if self.journal in ['MNRAS', 'A&A'] else None
            if 'ADS link' not in self.ads_scrapped_data[self.query]:
                self.ads_scrapped_data[self.query]['ADS link'] = join(self.base_links['ADS'], 'abs', self.query)
                self.ads_scrapped_data[self.query]['Publisher link via ADS gateway'] = gateway
                
            self.journal = self.query.split('.')[0][4:] if 'arXiv' not in self.query else 'arXiv'
            
            #Skip paper if not in commonly-used journals or paper is from arXIV
            if self.journal == 'arXiv' or not self.base:
                print('ACTION>>>>> Skipping since arXiv.\n')
                if 'Skipped' not in self.ads_scrapped_data[self.query]['Status']:
                    self.changed = True
                    self.ads_scrapped_data[self.query]['Status'] = 'Skipped'
                    self.ads_scrapped_data[self.query]['Action'] = 'Did not process since arXiv.' if self.journal == 'arXiv' else 'Did not process since only 1-2 papers in journal.'
                    self.save_overview()
                continue
                
            print('ADS gateway URL: %s\n'%gateway)
            #Request paper information from ADS gateway link to journal
            if self.rescan_online or 'Paper text' not in self.ads_scrapped_data[self.query] or self.ads_scrapped_data[self.query]['Status'] in ['Start', 'Failed']:
                self.changed = True
                response = self.sessions[self.base].get(gateway, headers = self.headers)
                sleep(self.slow_down_seconds_after_requests)
                got_access = self.check_passed(response)
                if got_access != 'Good':
                    self.ads_scrapped_data[self.query]['Status'] = 'Failed'
                    self.ads_scrapped_data[self.query]['Action'] = 'Aborted since: %s'%got_access
                else:
                    self.ads_scrapped_data[self.query]['Status'] = 'Page'
                    self.ads_scrapped_data[self.query]['Action'] = 'Loaded page'
                    self.ads_scrapped_data[self.query]['Paper text'] = response.text
                response.close()
            
            #Test access to journal
            if self.ads_scrapped_data[self.query]['Status'] in ['Failed', 'Skipped']:
                print('ACTION> %s: %s\n'%(self.ads_scrapped_data[self.query]['Status'], self.ads_scrapped_data[self.query]['Action']))
                self.save_overview()
                continue
            
            #Scrape paper data regarding tables
            if not self.rescan_online and 'Tables status' in self.ads_scrapped_data[self.query] and self.ads_scrapped_data[self.query]['Tables status'] in ['Scanned', 'Complete']: pass
            else:
                self.changed = True
                if self.journal in ['AJ', 'ApJ', 'ApJS', 'NAAS']:
                    try: table_links = self.scan_IOP_tables(self.ads_scrapped_data[self.query]['Paper text'])
                    except Exception as e:
                        print('IOP paper stopped since suspected has no figures.', e)
                        self.save_overview()
                        continue
                elif self.journal == 'MNRAS': table_links = self.scan_MNRAS_tables(self.ads_scrapped_data[self.query]['Paper text'])
                elif self.journal == 'A&A':
                    if self.base == 'A&A' and 'Tables at CDS' in self.ads_scrapped_data[self.query]['Paper text']:
                        '''self.ads_scrapped_data[self.query]['Status'] = 'Skipped'
                        self.ads_scrapped_data[self.query]['Action'] = 'Deferred to Vizier link'''
                        self.ads_scrapped_data[self.query]['Vizier link'] = BeautifulSoup(self.ads_scrapped_data[self.query]['Paper text']).find(href=True, title="Tables at CDS")['href']
                        print('NOTE> HAS Vizier link: %s\n'%self.ads_scrapped_data[self.query]['Vizier link'])
                        self.save_overview()
                        #continue
                    #else:
                    table_links = self.scan_AA_tables(self.ads_scrapped_data[self.query]['Paper text'])
                else: input('OOOPS...some journal slipped through the cracks. Check!!!')
                print('>>>>>Table links', table_links)
                self.populate_query_table(table_links)

            #Process table data through Pandas
            if self.ads_scrapped_data[self.query]['Tables status'] != 'Empty' and (self.rescan_online or self.redo_pandas or self.ads_scrapped_data[self.query]['Tables status'] in ['Scanned', 'Missing']):
                self.run_pandas()
                self.changed = True
                
            #Save overview information, including references to tables.
            if self.changed or self.rescan_online or self.redo_pandas: self.save_overview()
            
    def report_overal_stats(self):
        """Report overall success of scrape process"""
    
        status = {'Start':0, 'Incomplete':0, 'Page':0, 'Complete':0, 'Skipped':0, 'Failed':0}
        table_status = {'Empty':0, 'Missing':0, 'Scanned':0, 'Complete':0}
        for query in self.ads_scrapped_data.keys(): status[self.ads_scrapped_data[query]['Status']]+=1
        print('Overall Status>', status)
        for query in self.ads_scrapped_data.keys():
            try: table_status[self.ads_scrapped_data[query]['Tables status']]+=1
            except: print('No tables status, which is fine:', self.ads_scrapped_data[query]['Status'])
        print('Overall Status>', table_status)
            
    def save_overview(self):
        """Save overview query to file"""
        
        save_overview_path = join(self.base_directory, 'overviews', '%s-overview.txt'%self.query)
        with open(save_overview_path, 'w') as json_file: json.dump(self.ads_scrapped_data[self.query], json_file)
            
    def load_saved_overview(self):
        """Load overview query from file"""
        
        print(self.base_directory, 'overviews', '%s-overview.txt'%self.query)
        load_overview_path = join(self.base_directory, 'overviews', '%s-overview.txt'%self.query)
        if exists(load_overview_path):
            with open(load_overview_path, 'r') as json_file: self.ads_scrapped_data[self.query] = json.load(json_file)
        else: self.ads_scrapped_data[self.query] = {}
    
    def load_saved_tables(self):
        """Load table via overview keys"""
            
        self.ads_scrapped_tables[self.query] = {}
        if 'Table meta data' not in self.ads_scrapped_data[self.query]: pass
        else:
            for key in self.ads_scrapped_data[self.query]['Table meta data'].keys():
                if 'Pandas format' in self.ads_scrapped_data[self.query]['Table meta data'][key]:
                    for key2 in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'].keys():
                        if self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Status'] == 'Complete':
                            if key not in self.ads_scrapped_tables[self.query]: self.ads_scrapped_tables[self.query][key] = {}
                            if key2 not in self.ads_scrapped_tables[self.query][key]: self.ads_scrapped_tables[self.query][key][key2] = {}
                            self.ads_scrapped_tables[self.query][key][key2] = pd.read_json(join(self.base_directory, 'tables', '%s-set-%s-table-%s.text'%(self.query, key.replace('Table set ', ''), key2.replace('Table ',''))))
                
    def scan_AA_tables(self, response):
        """Scan an A&A page for tables"""
        
        return list(set([self.base_links[self.base] + link['href'] for link in BeautifulSoup(response).find_all(href=re.compile('T'), string=re.compile('Table *.'))]))
                
    def scan_IOP_tables(self, response):
        """Scan an IOPScience journal for tables"""
        
        #Parse table links
        table_links = list(set([self.base_links[self.base] + link['href'] for link in BeautifulSoup(response).find_all(href=re.compile('.txt'))]))
                
        #Older IOPScience links may use frames. Trying to step through a few other links to find tables if first pass fails.
        if len(table_links) > 0: pass
        else: table_links = self.scan_possible_IOP_outsource(response)
        
        return table_links
    
    def scan_possible_IOP_outsource(self, response):
        """Sometimes the journal points to an older page...which this function is designed to find some of the table formats within"""
        
        self.ads_scrapped_data[self.query]['IOP outsource trace'] = []
        #Traverse first link to frames
        traverse_link = self.base_links[self.base] + BeautifulSoup(response).find(href=re.compile('fulltext'))['href']
        print('IOP outsource linke>', traverse_link)
        
        #Get frames to link that links to table
        response_frames = self.sessions[self.base].get(traverse_link, headers = self.headers)
        got_access = self.check_passed(response_frames)
        self.ads_scrapped_data[self.query]['IOP outsource trace'].append({'Link': traverse_link, 'Response': response_frames.text, 'Status': 'Failed' if 'Good' not in access_status else 'Scanned', 'Action': ('Aborted since: %s'% access_status) if 'Good' not in access_status else 'Scanned link'})
        
        sleep(self.slow_down_seconds_after_requests)
        if got_access != 'Good':
            frame_article = join(traverse_link, BeautifulSoup(response_frames.text).find('frame', attrs={"name":"article"})['src'])
            response_frames.close()
        else:
            print('Failed to get IOP outsource frame link:', got_access)
            return []
        print('IOP frame link>', frame_article)
        
        #Traverse next frame to find table links
        deeper_frame = self.sessions[self.base].get(frame_article, headers = self.headers)
        got_access = self.check_passed(deeper_frame)
        self.ads_scrapped_data[self.query]['IOP outsource trace'].append({'Link': traverse_link, 'Response': deeper_frame.text, 'Status': 'Failed' if 'Good' not in access_status else 'Scanned', 'Action': ('Aborted since: %s'% access_status) if 'Good' not in access_status else 'Scanned link'})
        
        sleep(self.slow_down_seconds_after_requests)
        if got_access != 'Good':
            table_links = list(set([traverse_link + extension['href'] for extension in BeautifulSoup(deeper_frame.text).find_all('a', href=re.compile('.tb')) if extension['href'][-4:] == 'html']))
            deeper_frame.close()
            return table_links 
        else:
            print('Failed to get IOP outsource last link:', got_access)
            return []

    def scan_MNRAS_tables(self, response):
        """Scan MNRAS journal for tables"""
        
        #Parse table links
        return list(set([self.base_links[self.base] + extension['href'] for extension in  BeautifulSoup(response).find_all(href=re.compile('view-large'), class_="fig-view-orig openInAnotherWindow btn js-view-large")]))
        
    def populate_query_table(self, table_links):
        """Initialize table data dictionary and record table data prior to pandas processing. Will ignore if data already present."""
        
        #First set the paper text as the first place to scan for tables. NOTE repeats can happen between tables located in both paper text and contains an external link.
        if 'Table meta data' not in self.ads_scrapped_data[self.query]: self.ads_scrapped_data[self.query]['Table meta data'] = {}
        if 'Table set 0' not in self.ads_scrapped_data[self.query]['Table meta data']: self.ads_scrapped_data[self.query]['Table meta data']['Table set 0'] = {}
        self.ads_scrapped_data[self.query]['Table meta data']['Table set 0']['Link'] = self.ads_scrapped_data[self.query]['Publisher link via ADS gateway']
        self.ads_scrapped_data[self.query]['Table meta data']['Table set 0']['Response'] = self.ads_scrapped_data[self.query]['Paper text']
        self.ads_scrapped_data[self.query]['Table meta data']['Table set 0']['Status'] = 'Scanned' if self.ads_scrapped_data[self.query]['Status'] not in ['Failed', 'Skipped'] else 'Failed'
        self.ads_scrapped_data[self.query]['Table meta data']['Table set 0']['Action'] = 'Scanned link' if self.ads_scrapped_data[self.query]['Status'] not in ['Failed', 'Skipped'] else 'Aborted since main page not loaded'
        check = self.ads_scrapped_data[self.query]['Status'] in ['Failed', 'Skipped'] 
        
        #Loop through table links to add table data
        for index, table_link in enumerate(table_links):
            try:
                if len([key for key in self.ads_scrapped_data[self.query]['Table meta data'].keys() if key != 'Table set %s'%(index+1) and self.ads_scrapped_data[self.query]['Table meta data'][key]['Link'] == table_link]):
                    print('Duplicate link thus skipping')
                    continue
            except: print('Diction links not populated thus continuing to add link: ', table_link)
                
            if 'Table set %s'%(index+1) not in self.ads_scrapped_data[self.query]['Table meta data']: self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)] = {}
            if self.rescan_online or 'Status' not in self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)] or self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)]['Status'] == 'Failed':
                response = self.sessions[self.base].get(table_link, headers = self.headers)
                sleep(self.slow_down_seconds_after_requests)
                access_status = self.check_passed(response)
                
                self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)]['Link'] = table_link
                self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)]['Response'] = response.text
                self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)]['Status'] = 'Failed' if 'Good' not in access_status else 'Scanned'
                self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)]['Action'] = 'Scanned link' if 'Good' in access_status else ('Aborted since: %s'%access_status)
                if 'Failed' in self.ads_scrapped_data[self.query]['Table meta data']['Table set %s'%(index+1)]['Status']: check = True
            response.close()
        self.ads_scrapped_data[self.query]['Tables status'] = 'Empty' if len (table_links) == 0 else 'Missing' if check else 'Scanned'
                    
    def run_pandas(self):
        """Parse tables from web page via pandas"""
        
        if len(self.ads_scrapped_data[self.query]['Table meta data'].keys()) == 0: return
        
        done_ascii = []
        done_pandas = []
        table_unbroken = True
        for key in sorted(self.ads_scrapped_data[self.query]['Table meta data'].keys()):
            subtable_unbroken = True
            if self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] != 'Failed' and (self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] != 'Complete' or self.redo_pandas):
                print('>>>>>Table link', self.ads_scrapped_data[self.query]['Table meta data'][key]['Link'], self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'], '\n Query and key:', self.query, key)
                try:
                    #Determine type of table to be processed (ascii, HTML, ect)
                    response_text = self.ads_scrapped_data[self.query]['Table meta data'][key]['Response']
                    table_type = 'HTML' if '<table' in response_text else 'ASCII' if '-----' in response_text else 'local'
                    print('FIRST TRY on pandas', table_type)
                    if self.query in self.fix_messed_up_tables: response_text = response_text.replace(self.fix_messed_up_tables[self.query], '')
                    
                    #Run pandas to format table data
                    tables, extra = self.get_tables_using_html(response_text) if '<table' in response_text else self.get_tables_using_ascii(response_text) if '------' in response_text else self.get_tables_using_local(response_text)
                        
                    #Record possible table headers and footers
                    self.ads_scrapped_data[self.query]['Table meta data'][key]['Table captions or footers'] = np.array(extra).tolist()
                    
                    #Check if no tables exists
                    will_run = False
                    for pandas_table in tables:
                        will_run = True
                        break
                        
                    if not will_run:
                        table_unbroken = False
                        subtable_unbroken = False
                        print('EMPTY PANDAS TABLES', tables, response_text, self.ads_scrapped_data[self.query]['Table meta data'][key]['Link'], self.query, key)
                    
                    print('=====================================tables==================================', tables)
                    print('<<<<<<<<<EXTRA>>>>>>', extra)
                    
                    #Filter out duplicates
                    i2 = 0
                    duplicates = []
                    for pandas_table in tables:
                        #Check if duplicate
                        if table_type == 'ASCII':
                            duplicate = False
                            for pt in done_ascii:
                                try: duplicate = all(pt == pandas_table)
                                except:
                                    try: duplicate = pt == pandas_table
                                    except: duplicate = False
                                if duplicate: break
                            if duplicate: print('Skipping since table is a duplicate')
                            else:
                                done_ascii.append(pandas_table)
                                i2+=1
                        elif table_type != 'ASCII':
                            duplicate = False
                            for pt in done_pandas:
                                if pt == pandas_table.to_dict():
                                    duplicate = True
                                    break
                            if duplicate: print('Skipping since table is a duplicate')
                            else:
                                done_pandas.append(pandas_table.to_dict())
                                i2+=1
                        duplicates.append(duplicate)
                        if duplicate: continue

                        #Initialize table status data
                        if 'Pandas format' not in self.ads_scrapped_data[self.query]['Table meta data'][key]: self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'] = {}
                        if 'Table %s'%(i2+1) not in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']: self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)] = {'Status':'Start', 'Action':'Processing'}
                        if key not in self.ads_scrapped_tables[self.query]: self.ads_scrapped_tables[self.query][key] = {}
                        
                        if len(pandas_table) == 0:
                            table_unbroken = False
                            subtable_unbroken = False
                            self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)] = {'Status': 'Failed', 'Action': 'Aborted since empty or pandas failed to aquire data correctly.'}
                            print('Pandas table failed', pandas_table)
                            continue
                        elif self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)]['Status'] != 'Complete' or self.redo_pandas: self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)] = {'Status':'Start', 'Action':'Processing'}
                            
                        #Either read or load pandas table depanding on user specs.
                        try:
                            pandas_table_path = join(self.base_directory, 'tables', '%s-set-%s-table-%s.text'%(self.query, key.replace('Table set ',''), i2+1))
                            if not exists(pandas_table_path) or self.rescan_online or self.redo_pandas or self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)]['Status'] not in ['Failed', 'Complete']:
                                if table_type == 'ASCII': pandas_table.write(pandas_table_path, format = 'pandas.json')
                                else: pandas_table.to_json(pandas_table_path)
                                self.ads_scrapped_tables[self.query][key]['Table %s'%(i2+1)] = pandas_table
                                self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)] = {'Status': 'Complete', 'Action': 'Processed via pandas'}
                            elif self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)]['Status'] == 'Complete':
                                pandas_table_loaded = pd.read_json(pandas_table_path)
                                self.ads_scrapped_tables[self.query][key]['Table %s'%(i2+1)] = pandas_table_loaded
                            print('>Ran Pandas')
                        except Exception as e:
                            self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format']['Table %s'%(i2+1)] = {'Status': 'Failed', 'Action': 'Aborted since: %s'%e}
                            subtable_unbroken = False
                            table_unbroken = False
                            print('A pandas table failed', e, pandas_table)
                            
                    #Record if all tables were correctly processed
                    if subtable_unbroken:
                        if all(duplicates):
                            self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] = 'Skipped'
                            self.ads_scrapped_data[self.query]['Table meta data'][key]['Action'] = 'Not saving duplicates'
                        else:
                            self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] = 'Complete'
                            self.ads_scrapped_data[self.query]['Table meta data'][key]['Action'] = 'Scanned all subtables'
                    else:
                        self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] = 'Pandas failed'
                        self.ads_scrapped_data[self.query]['Table meta data'][key]['Action'] = 'At least a subtable failed pandas'

                except Exception as e:
                    print('Pandas for loop failed', e, self.ads_scrapped_data[self.query]['Table meta data'][key]['Response'])
                    subtable_unbroken = False
                    table_unbroken = False
                    self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] = 'Pandas failed'
                    self.ads_scrapped_data[self.query]['Table meta data'][key]['Action'] = 'Aborted since loop crash: %s'%e
                    input('on hold to check')
            elif self.ads_scrapped_data[self.query]['Table meta data'][key]['Status'] == 'Failed':
                table_unbroken = False
                print('>>>Broken link', self.ads_scrapped_data[self.query]['Table meta data'][key]['Link'], self.ads_scrapped_data[self.query]['Table meta data'][key]['Response'])
        #Record if all tables are correcly processed
        if table_unbroken:
            self.ads_scrapped_data[self.query]['Status'] = 'Complete'
            self.ads_scrapped_data[self.query]['Action'] = 'Scanned and processed all tables'
            self.ads_scrapped_data[self.query]['Tables status'] = 'Complete'
        else:
            self.ads_scrapped_data[self.query]['Status'] = 'Incomplete'
            self.ads_scrapped_data[self.query]['Action'] = 'One or more tables did not process'
            if self.ads_scrapped_data[self.query]['Tables status'] == 'Complete': self.ads_scrapped_data[self.query]['Tables status'] = 'Scanned'
                                    
    def check_passed(self, response):
        """Check if data loaded correctly"""
        
        try: 
            response.raise_for_status()
            if 'You do not currently have access to this article.' in response.text: return 'Dont have permission'
            if 'think that you are a bot' in response.text or 'confirm that you are not a robot' in response.text: return 'Blocked. Site thinks you are a bot. Possibly caused by overuse.'
        except Exception as e: return str(e)
        return 'Good'
        
    def get_tables_using_html(self, response):
        """Get table from HTML format"""
        
        return pd.read_html(response), {'Overview': caption.split('<p')[-1].split('</p')[0] for caption in response.split('caption') if '<table' in caption}
        
    def get_tables_using_ascii(self, response):
        """Load table via ascii"""
        
        return [Table.read(response, format = 'ascii')], {'Overview': response[:len(response) - response[::-1].index('---------------')]}
        
    def get_tables_using_local(self, response):
        """Get table using file write and load since pandas.load_html cannot read file via online"""
        
        #Set table and extra information containers
        tables = []
        extra_info = {}
        
        #Test if table is tabulated
        max_size = np.max([len(p.split('\t')) for p in response.split('\t\n')])
        if max_size >= 1:
            for index, trail_table in enumerate(response.split('\n\n')):
                #Record potential header and footer information vs table info depending on if \t is present in the row.
                if '\t' not in trail_table and trail_table != '': extra_info['Overview' if index == 0 else ('Info after table %s'%index)] = trail_table
                elif '\t' in trail_table:
                    max_size = np.max([len(p.split('\t')) for p in trail_table.split('\t\n')])
                    table_rows = ''
                    for suspect_table_row in trail_table.split('\t\n'):
                        if '\t' not in suspect_table_row:
                            if ('Overview data' if index == 0 else ('Info after table %s'%index)) not in extra_info: extra_info['Overview data' if index == 0 else ('Info after table %s'%index)] = suspect_table_row
                            else: extra_info['Overview data' if index == 0 else ('Info after table %s'%index)] += '\n%s'%suspect_table_row
                        else: table_rows += (suspect_table_row.replace('\n','') + '\t'*(max_size - len(suspect_table_row.split('\t'))) + '\n')
                    with open(join(self.base_directory, 'resources', 'temp.txt'), 'w') as file: file.write(table_rows)
                    tables.append(pd.read_table(join(self.base_directory, 'sandbox', 'temp.txt'), sep='\n', delimiter='\t', na_filter=False))
        else:
            print('Pandas local appeared to fail')
            tables.append(pd.read_table(join(self.base_directory, 'sandbox', 'temp.txt')))
            input('CHECK THIS')
            extra_info['Failed'] = 'Delimiter not anticipated'

        return tables, extra_info
        
    def get_standard_name_and_coords(self, table_row, map):
        """Set a coordiante based abbreviated summary name"""
        
        #Thankfully all papers appear to use J2000 format. Thus we do not check for J2000 here.
        if 'System Name' in map:
            if map['System Name'] == 'pdname': name = table_row.name
            else: name = table_row[map['System Name']]
        else: name = None
        
        standard_ra, standard_dec = self.convert_to_standard_ra_dec(table_row, map)
        if standard_ra is not None:
            #print('Computed Standard RA and DEC to derive a coordinate derived name')
            coord = SkyCoord(standard_ra, standard_dec, frame='fk5', unit='deg')
            ra_hour, ra_min, ra_sec = coord.ra.hms
            dec_sign, dec_degree, dec_arcmin, dec_arcsec = coord.dec.signed_dms
            dec_sign = '-' if dec_sign<0 else '+'
            ra_hms = "{ra_hour:02.0f}:{ra_min:02.0f}:{ra_sec:05.2f}".format(ra_hour=ra_hour,ra_min=ra_min,ra_sec=ra_sec)
            dec_dms = "{dec_sign}{dec_degree:02.0f}:{dec_arcmin:02.0f}:{dec_arcsec:05.2f}".format(dec_sign=dec_sign,dec_degree=dec_degree,dec_arcmin=dec_arcmin,dec_arcsec=dec_arcsec)
            standard_name = 'J%s%s' % (ra_hms.replace(':','')[:4], dec_dms.replace(':','')[:5])
            return standard_ra, standard_dec, standard_name
        elif name is not None and 'J' in name:
            if map['System Name'] == 'pdname': system_name = table_row.name
            else: system_name = table_row[map['System Name']]
            #print('Parsing coordinate name from system name (from table) that contains a "J" format:', system_name)
            system_name = 'J' + system_name.split('J')[1]
            
            #There are subtle character differences between signs. This code determines the sign used and splits the RA and DEC components from the system name.
            for sign in ['-', '–', '−', '+']:
                coords = system_name.split(sign)
                print('system trial', sign, coords)
                if len(coords) == 2: break
            ra, dec = coords
            ra = ra.replace('J','')[:4]
            dec = dec[:4]
            #Remove non-decimal related information
            non_decimal = re.compile(r'[^\d.]+')
            
            #Convert RA and DEC to a 4 digit form.
            ra = self.remove_non_numeric_related_formats(ra.split('.')[0])
            dec = self.remove_non_numeric_related_formats(dec.split('.')[0])
            print('Ra dec', ra, dec)
            system_name = 'J' + ra + sign + dec
            print('standard name', system_name)
            return '', '', system_name
        else: return '', '', ''
        
    def remove_non_numeric_related_formats(self, string, remove_plus = True):
        """Remove or replace numbers symbols that cause problems."""
    
        return "".join([n for n in string.replace('−','-').replace('+' if remove_plus else '','').replace('\xa0', '').replace('\u2009','').replace('...','') if n in '0123456789.+-: '])

    def convert_to_standard_ra_dec(self, table_row, map):
        """Convert input coordinate formats to a standard format"""
        
        #Use column map derived in inspection to convert table RA and DEC to a standard format
        if 'RA [°]' in map and table_row[map['RA [°]']] is not None:
            return float(self.remove_non_numeric_related_formats(str(table_row[map['RA [°]']]))), float(self.remove_non_numeric_related_formats(str(table_row[map['Dec [°]']])))
        elif 'RA in Hours' in map and table_row[map['RA in Hours']] is not None:
            RA_hour, RA_decimal = self.remove_non_numeric_related_formats(str(table_row[map['RA in Hours']])).split('.')
            dec_degree, dec_decimal = self.remove_non_numeric_related_formats(str(table_row[map['DEC in Degrees']])).split('.')
            return float(RA_hour)/24*360 + float(RA_decimal)/24*360, float(dec_degree) + float(dec_decimal)
        elif 'RA in Hours:Min:Sec' in map and table_row[map['RA in Hours:Min:Sec']] is not None:
            splitter = ':' if ':' in str(table_row[map['RA in Hours:Min:Sec']]) else ' '
            if splitter not in str(table_row[map['RA in Hours:Min:Sec']]):
                print('RA HMS test', self.parse_out_numbers(table_row[map['RA in Hours:Min:Sec']]))
                print('RA HMS result', self.parse_out_numbers(table_row[map['RA in Hours:Min:Sec']]))
                rcord = self.parse_out_numbers(table_row[map['RA in Hours:Min:Sec']])
                deccord = self.parse_out_numbers(table_row[map['Dec (+/-) Degree:Min:Sec']])
                if len(rcord) == 3: Rh, Rm, Rs = rcord
                elif len(rcord) == 4:
                    Rh, Rm, Rs, decimal = rcord
                    Rs += ('.' + decimal)
                else:
                    input('Failed RA coord')
                    return None, None
                if len(deccord) == 3: Dd, Dm, Ds = deccord
                elif len(deccord) == 4:
                    Dd, Dm, Ds, decimal = deccord
                    Ds += ('.' + decimal)
                else:
                    input('Failed DEC coord')
                    return None, None

                print('======splitter', splitter)
                print('Coords split', Rh, Rm, Rs, Dd, Dm, Ds)
            else:
                print('======splitter', splitter)
                Rh, Rm, Rs = self.remove_non_numeric_related_formats(str(table_row[map['RA in Hours:Min:Sec']])).split(splitter)
                Dd, Dm, Ds = self.remove_non_numeric_related_formats(str(table_row[map['Dec (+/-) Degree:Min:Sec']])).split(splitter)
                print(Rh, Rm, Rs, Dd, Dm, Ds)
            return float(Rh)/24*360 + float(Rm)/60 + float(Rs)/3600, float(Dd) + float(Dm)/60 + float(Ds)/3600
        elif 'RA in Degrees:Min:Sec' in map and table_row[map['RA in Degrees:Min:Sec']] is not None:
            splitter = ':' if ':' in str(table_row[map['RA in Degrees:Min:Sec']]) else ' '
            if splitter not in str(table_row[map['RA in Degrees:Min:Sec']]):
                Rh, Rm, Rs = self.parse_out_numbers(table_row[map['RA in Degrees:Min:Sec']])
                Dd, Dm, Ds = self.parse_out_numbers(table_row[map['Dec (+/-) Degree:Min:Sec']])
            else:
                Rd, Rm, Rs = self.remove_non_numeric_related_formats(str(table_row[map['RA in Degrees:Min:Sec']])).split(splitter)
                Dd, Dm, Ds = self.remove_non_numeric_related_formats(str(table_row[map['Dec (+/-) Degree:Min:Sec']])).split(splitter)
            return float(Rd) + float(Rm)/60 + float(Rs)/3600, float(Dd) + float(Dm)/60 + float(Ds)/3600
        elif 'RA (Hours part)' in map and table_row[map['RA (Hours part)']] is not None:
            R_degrees = self.remove_non_numeric_related_formats(str(table_row[map['RA (Hours part)']])) if 'RA (Hours part)' in map else self.remove_non_numeric_related_formats(str(table_row[map['RA (Degree part)']]))
            Rm, Rs = self.remove_non_numeric_related_formats(str(table_row[map['RA (Mins part)']])), self.remove_non_numeric_related_formats(str(table_row[map['RA (Secs part)']]))
            Dd, Dm, Ds = self.remove_non_numeric_related_formats(str(table_row[map['Dec (Degree part)']])), self.remove_non_numeric_related_formats(str(table_row[map['Dec (Arcmin part)']])), self.remove_non_numeric_related_formats(str(table_row[map['Dec (Arcsec part)']]))
            print(R_degrees, Rm, Rs, Dd, Dm, Ds)
            return float(R_degrees)/((24*360) if 'RA (Hours part)' in map else 1) + float(Rm)/60 + float(Rs)/3600, float(Dd) + float(Dm)/60 + float(Ds)/3600
        elif 'Position' in map and table_row[map['Position']] is not None:
            splitter = ':' if ':' in str(table_row[map['Position']]) else ' '
            for sign in ['-', '+', ' ']:
                coords = self.remove_non_numeric_related_formats(str(table_row[map['Position']]), remove_plus = False).split(sign)
                if len(coords) == 2: break
            Rh, Rm, Rs = coords[0].split(splitter)
            Dd, Dm, Ds = coords[1].split(splitter)
            Dd = sign + Dd
            return float(Rh)/24*360 + float(Rm)/60 + float(Rs)/3600, float(Dd) + float(Dm)/60 + float(Ds)/3600
        elif 'RA-Dec (Degrees)' in map and table_row[map['RA-Dec (Degrees)']] is not None:
            for sign in ['-', '+', ' ']:
                coords = self.remove_non_numeric_related_formats(str(table_row[map['RA-Dec (Degrees)']]), remove_plus = False).split(sign)
                if len(coords) == 2: break
            return float(self.remove_non_numeric_related_formats(coords[0])), float(sign + coords[1])
        elif 'MJD' in map and table_row[map['MJD']] is not None: return self.get_ra_dec_from_plate_mjd_fiberid(table_row[map['Plate']], table_row[map['MJD']], table_row[map['Fiberid']])
        elif 'MJD2' in map and table_row[map['MJD2']] is not None: return self.get_ra_dec_from_plate_mjd_fiberid(table_row[map['Plate2']], table_row[map['MJD2']], table_row[map['Fiberid2']])
        elif 'MJD3' in map and table_row[map['MJD3']] is not None: return self.get_ra_dec_from_plate_mjd_fiberid(table_row[map['Plate3']], table_row[map['MJD3']], table_row[map['Fiberid3']])
        else: return None, None
            
    def get_ra_dec_from_plate_mjd_fiberid(self, plate, mjd, fiberid):
        """Get RA and DEC via Plate, IFU, Fiberid"""
        
        #Look at dr16 Skyserver from SDSS for PLate, IFU, Fiberid info on RA and DEC
        sleep(1.5)
        session = Session()
        print('<<<Accessing from DR16>>>', r'https://skyserver.sdss.org/dr16/en/tools/search/x_results.aspx?searchtool=SQL&TaskName=Skyserver.Search.SQL&syntax=NoSyntax&ReturnHtml=true&cmd=--+This+query+does+a+table+JOIN+between+the+imaging+%28PhotoObj%29+and+spectra%0D%0A--+%28SpecObj%29+tables+and+includes+the+necessary+columns+in+the+SELECT+to+upload%0D%0A--+the+results+to+the+SAS+%28Science+Archive+Server%29+for+FITS+file+retrieval.%0D%0ASELECT+TOP+10%0D%0A+++s.ra%2Cs.dec%0D%0AFROM+SpecObj+AS+s%0D%0AWHERE+%0D%0A+++s.plate%3D' + str(plate) + r'+AND+s.' + r'mjd%3D' + str(mjd) + r'+AND+s.fiberid%3D' + str(fiberid) + r'%0D%0A%0D%0A&format=html&TableName=')

        sky_server_get_coords = pd.read_html(session.get(r'https://skyserver.sdss.org/dr16/en/tools/search/x_results.aspx?searchtool=SQL&TaskName=Skyserver.Search.SQL&syntax=NoSyntax&ReturnHtml=true&cmd=--+This+query+does+a+table+JOIN+between+the+imaging+%28PhotoObj%29+and+spectra%0D%0A--+%28SpecObj%29+tables+and+includes+the+necessary+columns+in+the+SELECT+to+upload%0D%0A--+the+results+to+the+SAS+%28Science+Archive+Server%29+for+FITS+file+retrieval.%0D%0ASELECT+TOP+10%0D%0A+++s.ra%2Cs.dec%0D%0AFROM+SpecObj+AS+s%0D%0AWHERE+%0D%0A+++s.plate%3D' + str(plate) + r'+AND+s.' + r'mjd%3D' + str(mjd) + r'+AND+s.fiberid%3D' + str(fiberid) + r'%0D%0A%0D%0A&format=html&TableName=', headers = self.headers).text)[0].iloc[1]
        session.close()
        return float(sky_server_get_coords[0]), float(sky_server_get_coords[1])
                    
    def set_system_data(self):
        """Set data per system entry in a table"""
    
        for index, self.query in enumerate(self.bibcodes[self.end:self.start:-1]):
            print('INDEX>>>', index)
            
            referenced = self.check_referance_added_to_MLD_references(self.query, self.ads_scrapped_data[self.query]['Paper Overview'])
            if referenced:
                print('Paper already cited: %s. Skipping'%self.query)
                continue
            else: pass
            
            #Get data from each table via the table column map determined in table inspection
            if ('Status' in self.ads_scrapped_data[self.query] and self.ads_scrapped_data[self.query]['Status'] == 'Complete'):
                for key in sorted(self.ads_scrapped_data[self.query]['Table meta data'].keys()):
                    if 'Pandas format' in self.ads_scrapped_data[self.query]['Table meta data'][key]:
                        print('LINK>>>>', self.ads_scrapped_data[self.query]['Table meta data'][key]['Link'])
                        for key2 in sorted(self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'].keys()):
                            if 'Inspection' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2] and 'Table map to MasterLens database' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']:
                                for row_index in range(0 if 'Append first table row to header' not in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'] else 1, self.ads_scrapped_tables[self.query][key][key2].shape[1 if 'Table columns are horizontal' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'] else 0]):
                                    #Get data per table row depending on if table is rotated or not
                                    if 'Table columns are horizontal' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']: table_row = self.ads_scrapped_tables[self.query][key][key2].T.iloc[row_index]
                                    else:
                                        table_row = self.ads_scrapped_tables[self.query][key][key2].iloc[row_index]
                                    map = self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Table map to MasterLens database']
                                                                        
                                    if 'Condition' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Table map to MasterLens database'] and 'Word to recognize candidate in column' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'] and self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']['Word to recognize candidate in column']:
                                        filter_string = self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']['Word to recognize candidate in column']
                                        description = table_row.name if map['Condition'] == 'pdname' else table_row[map['Condition']]
                                        if filter_string.startswith('!'):
                                            filter_string = filter_string.split('!')[-1]
                                            if filter_string in description:
                                                print('Skipping row since condition not met', table_row)
                                                continue
                                        elif filter_string not in description:
                                            print('Skipping row since condition not met', table_row)
                                            continue
                                            
                                    print('Table', self.ads_scrapped_tables[self.query][key][key2], map)
                                    
                                    
                                    #Some tables have repeats forms horizontally. This condition scans horizontal repeats.
                                    if 'Repeats format horizontally' in self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']:
                                        sub_map = {}
                                        splitter = int(self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes']['Repeats format horizontally'])
                                        for i in range(0, len(table_row), splitter):
                                            #Map horizontal repeated form
                                            for table_key in table_row[i:i+splitter].keys():
                                                print('table_key', table_key)
                                                for map_key in map:
                                                    print('map_key', map_key, map[map_key], table_key)
                                                    if map[map_key] == table_key: sub_map[map_key] = table_key
                                            self.set_object_data(key, key2, table_row, sub_map, self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'])
                                    else: self.set_object_data(key, key2, table_row, map, self.ads_scrapped_data[self.query]['Table meta data'][key]['Pandas format'][key2]['Inspection']['Notes'])
                                    
    def set_object_data(self, key, key2, table_row, map, action_map):
        """Get system data from a table"""
        
       
        #These keys cause integer vs string issues in dictionary reference translations.
        oversimplified_keys = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']
        
        #I have to converted several tables ids to strings during inspection process, which I reverse here. Can improve when have time to fix.
        print('map===>',map)
        print('table_row====>', table_row)
        print('table_row keys====>', table_row.keys())
        print('action_map====>', action_map)
        empty = True
        for mkey in map:
            if map[mkey] in oversimplified_keys: new_value = int(map[mkey])
            else: new_value = map[mkey]
            map[mkey] = new_value
            print('>>>>', mkey, map[mkey], new_value)
            try:
                if new_value == 'pdname' and table_row.name is not None: empty = False
                elif table_row[new_value.split('@')[0] if '@' in str(new_value) else new_value] is not None: empty = False
            except:
                if table_row[int(new_value.split('@')[0]) if '@' in str(new_value) else int(new_value)] is not None: empty = False
        if empty: print('Skipping row since empty')
        else:
            #Convert to a standardized format
            
            try: standard_ra, standard_dec, standard_name = self.get_standard_name_and_coords(table_row, map)
            except Exception as e:
                print(self.ads_scrapped_tables[self.query][key][key2], '\n>Standardize system Failed:', e)
                standard_name, standard_ra, standard_dec = '', '', ''
                print('Problem with data:', table_row, map)
                #testi = input('Retry to see bug? (type y for yes):')
                #if testi == 'y': standard_ra, standard_dec, standard_name = self.get_standard_name_and_coords(table_row, map)


            if 'Cluster Sources Table' in action_map and standard_ra != '':
                if 'Word to recognize name is of lens and NOT source' in action_map:
                    if action_map['Word to recognize name is of lens and NOT source'] in table_row[map['Source names']]:
                        self.cluster_lens_name = table_row[map['Source names']]
                        standard_name = self.cluster_lens_name + ''
                    else: standard_name = self.cluster_lens_name + '[' + table_row[map['Source names']] + ']'
                elif 'Name,Ra,Dec of cluster or group lens' in action_map:
                    self.cluster_lens_name = action_map['Name,Ra,Dec of cluster or group lens']
                    standard_name = self.cluster_lens_name + '[' + (str(table_row[map['Source names']]) if 'Source names' in map else standard_name) + ']'


            if standard_name is '': print('Could not define system', table_row, map)
            else:
                if standard_name not in self.lens_objects: self.lens_objects[standard_name] = {}
                if 'Cluster Sources Table' in action_map:
                    if 'System Name' in self.lens_objects[standard_name]: self.lens_objects[standard_name]['System Name'].append({'value': standard_name, 'tracer': {'update status': 'is or in cluster', 'weight':5}})
                    else: self.lens_objects[standard_name]['System Name'] = [{'value': standard_name, 'tracer': {'update status': 'is or in cluster', 'weight':5}}]
                if 'Standard RA' not in self.lens_objects: self.lens_objects[standard_name]['Standard RA'] = []
                if 'Standard DEC' not in self.lens_objects: self.lens_objects[standard_name]['Standard DEC'] = []

                #Save standard coordinate info
                self.lens_objects[standard_name]['Standard RA'].append({'value': standard_ra, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': 'Not yet included', 'weight':0}})
                self.lens_objects[standard_name]['Standard DEC'].append({'value': standard_dec, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': 'Not yet included', 'weight':0}})

                self.set_coord_details(standard_name, 0, key, key2, 'Not yet included')
                #Save reference information via a conversion from ADS bibform to MLD bibform
                if 'References' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['References'] = [self.ads_to_mld_reference_interpreter[self.query]]
                elif self.ads_to_mld_reference_interpreter[self.query] not in self.lens_objects[standard_name]['References']: self.lens_objects[standard_name]['References'].append(self.ads_to_mld_reference_interpreter[self.query])

                if 'System Name' in self.lens_objects[standard_name]:
                    print(self.lens_objects[standard_name]['System Name'])
                    system_names = [entry['value'] for entry in self.lens_objects[standard_name]['System Name']]
                else: system_names = []
                system_names.append(standard_name)

                self.lens_objects[standard_name]['References'] = self.get_all_papers_referenced(self.lens_objects[standard_name]['References'], system_names)

                #I wish this could work but there is no standard format that is robust across tables. Thus this method is deprecated.
                '''if 'References' in map:
                    for reference in table_row[map['References']]:
                        if reference not in self.lens_objects[standard_name]['References']: self.lens_objects[standard_name]['References'].append(reference)'''

                #Save any related method and error specified from inspection map of columns
                for mkey in map.keys():
                    if not isinstance(map[mkey], int) and '@' in map[mkey]:
                        mvalue, method, merror = map[mkey].split('@')
                        if mvalue in oversimplified_keys: mvalue = int(mvalue)
                        if merror in oversimplified_keys: merror = int(merror)
                        if mvalue == merror:
                            entry = self.parse_out_numbers(table_row[mvalue])
                            if len(entry) == 2:
                                value, error = entry
                                print(mvalue, entry, table_row[mvalue], table_row, map)
                            elif len(entry) == 3:
                                value = entry[0]
                                error = str((abs(float(entry[1])) + abs(float(entry[2])))/2)  
                                print(mvalue, entry, table_row[mvalue], table_row, map)
                            else:
                                print(mvalue, entry, table_row[mvalue], table_row, map)
                                input('Paused for you to check error parsing')
                                try:
                                    value = str(float(entry[0]))
                                    error = ''
                                except:
                                    value = ''
                                    error = ''
                        else:
                            value = table_row[mvalue]
                            error = table_row[merror] if merror != '' else ''
                    elif map[mkey] == 'pdname': value, method, error = table_row.name, '', ''
                    else: value, method, error = table_row[map[mkey]], '', ''

                    #Attempt a weight system used to determine which value is best to use for MLD. Check how we want to improve this!!!
                    weight = 0
                    weight += (1 if error != '' else 0)
                    weight += (1 if method in ["spectroscopic", "SIE model", "MLD"] else -1 if method in ["Reference redshift", ""] else 0)

                    #Save system information. The "replace(2,'') is to allow inspection of repeated column types (i.e. Plate, Plate2, Plate3)
                    if value not in ['', 'NaN', ' NaN', None]:
                        if mkey not in self.lens_objects[standard_name]: self.lens_objects[standard_name][mkey.replace('2','').replace('3','')] = []
                        self.lens_objects[standard_name][mkey.replace('2','').replace('3','')].append({'value': str(value).replace(' NaN','').replace('NaN',''), 'method': str(method).replace(' NaN','').replace('NaN',''), 'error': str(error).replace(' NaN','').replace('NaN',''), 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': 'Not yet included', 'weight':weight}})
                        print('lens object values', self.lens_objects[standard_name])
                if 'Lens type' in action_map: self.lens_objects[standard_name]['Lens type'] = action_map['Lens type']
                if 'Discovery' in action_map: self.lens_objects[standard_name]['Discovery'] = action_map['Discovery']
                 
    def set_coord_details(self, standard_name, tracer, weight, key, key2, update_status):
        if self.lens_objects[standard_name]['Standard RA']:
            coord = SkyCoord(self.lens_objects[standard_name]['Standard RA'], self.lens_objects[standard_name]['Standard Dec'], frame='fk5', unit='deg')
            
            hour, mn, sec = coord.ra.hms
            if 'RA (Hours part)' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['RA (Hours part)'] = []
            if 'RA (Mins part)' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['RA (Mins part)'] = []
            if 'RA (Secs part)' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['RA (Secs part)'] = []
            self.lens_objects[standard_name]['RA (Hours part)'].append({'value': hour, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
            self.lens_objects[standard_name]['RA (Mins part)'].append({'value': mn, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
            self.lens_objects[standard_name]['RA (Secs part)'].append({'value': sec, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
            
            
            if 'Dec (Degree part)' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['Dec (Degree part)'] = []
            if 'Dec (Arcmin part)' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['Dec (Arcmin part)'] = []
            if 'Dec (Arcsec part)' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['Dec (Arcsec part)'] = []
            dec_sign,deg, mn, sec = coord.dec.signed_dms
            dec_sign = '-' if dec_sign<0 else '+'
            dec_sign, self.lens_objects[standard_name]['Dec (Degree part)'].append({'value': dec_sign+deg, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
            self.lens_objects[standard_name]['Dec (Arcmin part)'].append({'value': mn, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
            self.lens_objects[standard_name]['Dec (Arcsec part)'].append({'value': sec, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
            
            if 'Dec [°]' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['Dec [°]'] = []
            if 'RA [°]' not in self.lens_objects[standard_name]: self.lens_objects[standard_name]['RA [°]'] = []
            self.lens_objects[standard_name]['Dec [°]'].append({'value': coord.dec.deg, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
            self.lens_objects[standard_name]['RA [°]'].append({'value': coord.ra.deg, 'tracer': {'bibcode':self.ads_to_mld_reference_interpreter[self.query], 'table set': key, 'table': key2, 'update status': update_status, 'weight':weight}})
                        
    def write_pdfs(self):
        """Write pdfs of each paper"""
        
        for self.query in self.bibcodes:            
            self.load_saved_overview()
            self.load_saved_tables()
            pdf_path = join(self.base_directory, 'pdfs', '%s.pdf'%self.query)
            
            if 'Publisher link via ADS gateway' in self.ads_scrapped_data[self.query] and not exists(pdf_path):
                session = Session()
                response = session.get(self.ads_scrapped_data[self.query]['Publisher link via ADS gateway'].replace('PUB_HTML', 'PUB_PDF'), headers = self.headers)
                check_response = self.check_passed(response)
                if 'Good' in check_response:
                    with open(pdf_path, 'wb') as file: file.write(response.content)
                    print('Done:', self.query)
                else: print('Failed to get pdf for', self.query)
                session.close()
                del response
                
    def load_MLD_references(self):
        """load Masterlens database reference entries from xml file."""
        
        #Set keys including a few manual added forms that have potentially bad bibcode references.
        self.ads_to_mld_reference_interpreter = {self.convert_to_mld_reference_form('1993LIACo..31..153S'): 'LIAC31(1993)', self.convert_to_mld_reference_form('2007NJPh....9..447J'):'New Journal of Physics9(2007)447', self.convert_to_mld_reference_form('2008MNRAS.386.2065C'): 'FORJ0332-3557-2', 'IDarkMatter(2008)':'IDarkMatter(2008)', self.convert_to_mld_reference_form('2007astro.ph..1537S'):'astro-ph/0701537', self.convert_to_mld_reference_form('2019ApJS..243...14O'):'http://arxiv.org/abs/1904.02806', self.convert_to_mld_reference_form('2002astro.ph..4513F'):'astro-ph/0204513', self.convert_to_mld_reference_form('2001astro.ph..2340K'):'astro-ph/0102340', self.convert_to_mld_reference_form('2001astro.ph..2341K'):'astro-ph/0102341', 'MNRAS472(2017)4038A':'arXiv:1702.00406', 'arXiv1110(2011)3784': '2011arXiv:1110.3784', 'ApJ774(2013)124': 'arXiv:1306.5240', 'ApJ764(2013)186': 'arXiv:1207.5776', 'MNRAS429(2013)482': 'arXiv:1209.0458', 'ApJ758(2012)26': 'arXiv:1202.1645', 'MNRAS439(2014)3392': 'arXiv:1206.3412', 'ApJ765(2013)139': 'arXiv:1206.2011', 'ApJ777(2013)97': 'arXiv1307.4764', 'ApJ777(2013)98': 'arXiv1307.4759'}
        
        self.ads_to_mld_reference_interpreter_inverted = {self.ads_to_mld_reference_interpreter[key]:key for key in self.ads_to_mld_reference_interpreter}
        
        #Set containers used to determine if a reference needs to be uploaded to MLD and what is the reference information.
        self.done_references_bib = [key for key in self.ads_to_mld_reference_interpreter.keys()]
        self.done_references_title = []
        self.mld_reference = {}
        self.reference_id = {}
        self.update_reference = []
        
        #Process xml reference file.
        doc = xml.dom.minidom.parse(join(self.base_directory, 'resources', 'references.xml'))
        mld_references = doc.getElementsByTagName('reference')
        
        for mld_reference in mld_references[:]:
            #First standardize the bibcode form
            print('PRE citation', mld_reference.getElementsByTagName('identifier')[0].childNodes[0].nodeValue)
            try: url = mld_reference.getElementsByTagName('adsabs_link')[0].childNodes[0].nodeValue
            except: url = None
            if mld_reference.getElementsByTagName('identifier')[0].childNodes[0].nodeValue in self.ads_to_mld_reference_interpreter_inverted: citation = self.ads_to_mld_reference_interpreter_inverted[mld_reference.getElementsByTagName('identifier')[0].childNodes[0].nodeValue]
            else: citation = self.convert_to_mld_reference_form(mld_reference.getElementsByTagName('identifier')[0].childNodes[0].nodeValue, url)
            print('>>>Bibcode from MLD:',citation)
            self.ads_to_mld_reference_interpreter[citation] = mld_reference.getElementsByTagName('identifier')[0].childNodes[0].nodeValue
            if citation not in self.done_references_bib: self.done_references_bib.append(citation)
            if mld_reference.getElementsByTagName('title')[0].childNodes[0].nodeValue not in self.done_references_title: self.done_references_title.append(mld_reference.getElementsByTagName('title')[0].childNodes[0].nodeValue)
            
            #Now record the bib information
            query = mld_reference.getElementsByTagName('identifier')[0].childNodes[0].nodeValue
            self.mld_reference[query] = {}
            self.mld_reference[query] = {"action":"Do not upload", "identifier":query, "author": mld_reference.getElementsByTagName('author')[0].childNodes[0].nodeValue, "title": mld_reference.getElementsByTagName('title')[0].childNodes[0].nodeValue, "journalID": mld_reference.getAttribute('referenceID')}
            
            try: self.mld_reference[query]["journal"] = mld_reference.getElementsByTagName('journal')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["abstract"] = mld_reference.getElementsByTagName('abstract')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["number"] = mld_reference.getElementsByTagName('number')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["year"] = mld_reference.getElementsByTagName('year')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["month"] = mld_reference.getElementsByTagName('month')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["pages"] = mld_reference.getElementsByTagName('pages')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["volume"] = mld_reference.getElementsByTagName('volume')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["doi"] = mld_reference.getElementsByTagName('doi_link')[0].childNodes[0].nodeValue
            except: pass
            try: self.mld_reference[query]["ads"] = mld_reference.getElementsByTagName('adsabs_link')[0].childNodes[0].nodeValue
            except: pass
            self.reference_id[query] = mld_reference.getAttribute('referenceID')
                
                
    def load_MLD_kinds_ids(self):
        """Load Masterlens Database lens kinds via kinds.xml"""
        
        doc = xml.dom.minidom.parse(join(self.base_directory, 'resources', 'kinds.xml'))
        kinds = doc.getElementsByTagName('kind')
        
        self.lens_type_id = {}
        for kind in kinds[:]:
            try: self.lens_type_id[kind.getElementsByTagName('acronym')[0].childNodes[0].nodeValue] = kind.getAttribute('kindID')
            except: pass #This try statement is here to avoid repeated nested kind statements, with the outer being empty.
        self.lens_type_id_inverted = {self.lens_type_id[key]:key for key in self.lens_type_id}
        
        
    def load_MLD_surveys_ids(self):
        """Load Masterlens Database discovery surveys via discovery.xml"""
        
        doc = xml.dom.minidom.parse(join(self.base_directory, 'resources', 'discoveries.xml'))
        kinds = doc.getElementsByTagName('discovery')
        
        self.discovery_id = {}
        for discovery in kinds[:]: self.discovery_id[discovery.getElementsByTagName('acronym')[0].childNodes[0].nodeValue] = discovery.getAttribute('discoveryID')
        self.discovery_id_inverted = {self.discovery_id[key]:key for key in self.discovery_id}
                                  
                
    def load_MLD_lenses(self):
        """Load Masterlens Database lenses via lenses.xml"""
        
        doc = xml.dom.minidom.parse(join(self.base_directory, 'resources', 'lenses.xml'))
        lenses = doc.getElementsByTagName('lens')
        
        self.lens_objects = {}
        for mld_entry in lenses[:]:
            standard_ra, standard_dec, standard_name = self.get_standard_name_and_coords({'RA [°]': str(mld_entry.getElementsByTagName('ra_coord')[0].childNodes[0].nodeValue), 'Dec [°]': str(mld_entry.getElementsByTagName('dec_coord')[0].childNodes[0].nodeValue)}, {'RA [°]':'RA [°]', 'Dec [°]':'Dec [°]'})
            if standard_name not in self.lens_objects: self.lens_objects[standard_name] = {'System Name':[], 'Discovery Date':[], 'RA (Hours part)':[], 'RA (Mins part)':[], 'RA (Secs part)':[], 'RA [°]': [], 'Dec (Degree part)': [], 'Dec (Arcmin part)': [], 'Dec (Arcsec part)': [], 'Dec [°]': [], 'Lens Grade': [], 'Number of images': [], 'Einstein_R ["]': [], 'z_Lens': [], 'z_Lens error':[], 'z_Lens quality':[], 'z_Source error':[], 'z_Source quality':[], 'z_Source(s)': [], 'Stellar velocity disp': [], 'Standard RA':[], 'Standard DEC':[], 'MLD_ID':[], 'Description':[], 'Lens type':[], 'Lens type MLD_ID':[], 'Discovery':[], 'Discovery_MLD_ID':[], 'MLD SDSS link':[], 'MLD ADSABS link':[], 'MLD NED link':[], 'MLD APOD link':[], 'References_MLD_ID':[], 'Has external link for SDSS':[], 'Has external link for ADSABS':[], 'Has external link for NED':[], 'Has external link for APOD':[]}
            
            self.lens_objects[standard_name]['MLD_ID'].append({'value': mld_entry.getAttribute('lensID'), 'tracer': {'update status': 'in MLD', 'weight':10}})
            self.lens_objects[standard_name]['Lens type MLD_ID'].append({'value': mld_entry.getElementsByTagName('kind')[0].getAttribute('kindID'), 'tracer': {'update status': 'in MLD', 'weight':10}})
            self.lens_objects[standard_name]['Discovery_MLD_ID'].append({'value': mld_entry.getElementsByTagName('discovery')[0].getAttribute('discoveryID'), 'tracer': {'update status': 'in MLD', 'weight':10}})
            self.lens_objects[standard_name]['Standard RA'].append({'value': standard_ra, 'tracer': {'update status': 'in MLD', 'weight':10}})
            self.lens_objects[standard_name]['Standard DEC'].append({'value': standard_dec, 'tracer': {'update status': 'in MLD', 'weight':10}})
            self.set_coord_details(standard_name, 10, 'MLD', 'MLD', 'in MLD')
            mld_entry_load_types = {'System Name': 'system_name', 'Discovery Date': 'discovery_date', 'Description': 'description', 'Lens type': 'kind', 'Discovery': 'discovery', 'RA (Hours part)': 'ra_hrs', 'RA (Mins part)': 'ra_mins', 'RA (Secs part)': 'ra_secs', 'RA [°]': 'ra_coord', 'Dec (Degree part)': 'dec_degrees', 'Dec (Arcmin part)': 'dec_arcsec', 'Dec (Arcsec part)': 'dec_arcsec', 'Dec [°]': 'dec_coord', 'Lens Grade': 'lensgrade', 'Number of images': 'number_images'}
            for tkey in mld_entry_load_types.keys():
                try: self.lens_objects[standard_name][tkey].append({'value': mld_entry.getElementsByTagName(mld_entry_load_types[tkey])[0].childNodes[0].nodeValue, 'tracer': {'update status': 'in MLD', 'weight':1 if tkey in ['Lens Grade', 'Number of images'] else 10}})
                except: pass
            
            #Check if external links are present
            try:
                self.lens_objects[standard_name]['MLD SDSS link'].append({'value': mld_entry.getElementsByTagName('sdss_link')[0].childNodes[0].nodeValue, 'tracer': {'update status': 'in MLD', 'weight':1}})
                self.lens_objects[standard_name]['Has external link for SDSS'].append({'value': '1', 'tracer': {'update status': 'in MLD', 'weight':1}})
            except: pass
            try:
                self.lens_objects[standard_name]['MLD ADSABS link'].append({'value': mld_entry.getElementsByTagName('adsabs_link')[0].childNodes[0].nodeValue, 'tracer': {'update status': 'in MLD', 'weight':10}})
                self.lens_objects[standard_name]['Has external link for ADSABS'].append({'value': '1', 'tracer': {'update status': 'in MLD', 'weight':10}})
            except: pass
            try:
                self.lens_objects[standard_name]['MLD NED link'].append({'value': mld_entry.getElementsByTagName('ned_link')[0].childNodes[0].nodeValue, 'tracer': {'update status': 'in MLD', 'weight':1}})
                self.lens_objects[standard_name]['Has external link for NED'].append({'value': '1', 'tracer': {'update status': 'in MLD', 'weight':10}})
            except: pass
            try:
                self.lens_objects[standard_name]['MLD APOD link'].append({'value': mld_entry.getElementsByTagName('apod_link')[0].childNodes[0].nodeValue, 'tracer': {'update status': 'in MLD', 'weight':1}})
                self.lens_objects[standard_name]['Has external link for APOD'].append({'value': '1', 'tracer': {'update status': 'in MLD', 'weight':10}})
            except: pass
            
            #Record values with potential errors
            values_with_errors_dict = {'theta_e':'Einstein_R ["]', 'z_lens':'z_Lens', 'z_source': 'z_Source(s)', 'vdisp': 'Stellar velocity disp'}
            for key in values_with_errors_dict.keys():
                weight = 0
                try:
                    value = mld_entry.getElementsByTagName(key)[0].childNodes[0].nodeValue
                except: value = ''
            
                try:
                    error = mld_entry.getElementsByTagName(key + '_err')[0].childNodes[0].nodeValue
                    weight = 2 if float(error) < (.005 if key in ['z_lens', 'z_source'] else .1 if key == 'theta_e' else 50) else 1
                except:
                    error = ''
                    weight = -1
                
                self.lens_objects[standard_name][values_with_errors_dict[key]].append({'value': value, 'method':'MLD', 'error': error, 'tracer': {'update status': 'in MLD', 'weight':weight}})
            self.lens_objects[standard_name]['References'] = [r.firstChild.data for r in mld_entry.getElementsByTagName('reference')]
            self.lens_objects[standard_name]['References_MLD_ID'] = [r.getAttribute('referenceID') for r in mld_entry.getElementsByTagName('reference')]
        
    def load_sugohi(self):
        sugohi_key={'SuGOHI1':'PASJ70S(2018)29S', 'SuGOHI2':'ApJ867(2018)107W', 'SuGOHI3':'A&A630A(2019)71S', 'SuGOHI4':'A&A636A(2020)87C', 'SuGOHI5':'MNRAS495(2020)1291J', 'SuGOHI6':'A&A642A(2020)148S', 'SuGOHI7':'MNRAS502(2021)1487J'}
        with open('list_public.csv', newline='') as csvfile:
            sugohi = csv.reader(csvfile, delimiter=',')
            for row in sugohi:
                if 'Name' in row: continue
                standard_ra, standard_dec, standard_name = self.get_standard_name_and_coords({'RA [°]': row[1], 'Dec [°]': row[2]}, {'RA [°]':'RA [°]', 'Dec [°]':'Dec [°]'})
                if standard_name not in self.lens_objects: self.lens_objects[standard_name] = {'System Name':[], 'Discovery Date':[], 'RA (Hours part)':[], 'RA (Mins part)':[], 'RA (Secs part)':[], 'RA [°]': [], 'Dec (Degree part)': [], 'Dec (Arcmin part)': [], 'Dec (Arcsec part)': [], 'Dec [°]': [], 'Lens Grade': [], 'Number of images': [], 'Einstein_R ["]': [], 'z_Lens': [], 'z_Lens error':[], 'z_Lens quality':[], 'z_Source error':[], 'z_Source quality':[], 'z_Source(s)': [], 'Stellar velocity disp': [], 'Standard RA':[], 'Standard DEC':[], 'MLD_ID':[], 'Description':[], 'Lens type':[], 'Lens type MLD_ID':[], 'Discovery':[], 'Discovery_MLD_ID':[], 'MLD SDSS link':[], 'MLD ADSABS link':[], 'MLD NED link':[], 'MLD APOD link':[], 'References_MLD_ID':[], 'Has external link for SDSS':[], 'Has external link for ADSABS':[], 'Has external link for NED':[], 'Has external link for APOD':[]}
                self.lens_objects[standard_name]['Standard RA'].append({'value': standard_ra, 'tracer': {'update status': 'in SuGOHI', 'weight':5}})
                self.lens_objects[standard_name]['Standard DEC'].append({'value': standard_dec, 'tracer': {'update status': 'in SuGOHI', 'weight':5}})
                self.set_coord_details(standard_name, 5, 'SuGOHI', 'SuGOHI', 'in SuGOHI')
                self.lens_objects[standard_name]['z_Lens'].append({'value': row[3], 'tracer': {'update status': 'in SuGOHI', 'weight':5}})
                self.lens_objects[standard_name]['z_Source(s)'].append({'value': row[4], 'tracer': {'update status': 'in SuGOHI', 'weight':5}})
                self.lens_objects[standard_name]['System Name'].append({'value': 'HSC'+standard_name, 'tracer': {'update status': 'in SuGOHI', 'weight':5}})
                for sr in row[-1].split(' '):
                    if 'SuGOHI' in sr: break
                if 'References' in self.lens_objects[standard_name]: self.lens_objects[standard_name]['References'].append(sr)
                else: self.lens_objects[standard_name]['References'] = [sr]
                
            
           
    def load_silo_eboss(self):
        hdu = fits.open('silo_eboss_detections-1.0.1.fits')
        for candidate in hdu['DETECTIONS'].data:
            standard_ra, standard_dec, standard_name = self.get_standard_name_and_coords({'RA [°]': candidate['RA'], 'Dec [°]': candidate['DEC']}, {'RA [°]':'RA [°]', 'Dec [°]':'Dec [°]'})
            if standard_name not in self.lens_objects: self.lens_objects[standard_name] = {'System Name':[], 'Discovery Date':[], 'RA (Hours part)':[], 'RA (Mins part)':[], 'RA (Secs part)':[], 'RA [°]': [], 'Dec (Degree part)': [], 'Dec (Arcmin part)': [], 'Dec (Arcsec part)': [], 'Dec [°]': [], 'Lens Grade': [], 'Number of images': [], 'Einstein_R ["]': [], 'z_Lens': [], 'z_Lens error':[], 'z_Lens quality':[], 'z_Source error':[], 'z_Source quality':[], 'z_Source(s)': [], 'Stellar velocity disp': [], 'Standard RA':[], 'Standard DEC':[], 'MLD_ID':[], 'Description':[], 'Lens type':[], 'Lens type MLD_ID':[], 'Discovery':[], 'Discovery_MLD_ID':[], 'MLD SDSS link':[], 'MLD ADSABS link':[], 'MLD NED link':[], 'MLD APOD link':[], 'References_MLD_ID':[], 'Has external link for SDSS':[], 'Has external link for ADSABS':[], 'Has external link for NED':[], 'Has external link for APOD':[]}
            self.lens_objects[standard_name]['Standard RA'].append({'value': standard_ra, 'tracer': {'update status': 'in SILO', 'weight':10}})
            self.lens_objects[standard_name]['Standard DEC'].append({'value': standard_dec, 'tracer': {'update status': 'in SILO', 'weight':10}})
            self.set_coord_details(standard_name, 5, 'SILO', 'SILO', 'in SILO')
            self.lens_objects[standard_name]['z_Lens'].append({'value': candidate['Z_NOQSO'], 'tracer': {'update status': 'in SILO', 'weight':5}})
            self.lens_objects[standard_name]['z_Lens error'].append({'value': candidate['ZERR_NOQSO'], 'tracer': {'update status': 'in SILO', 'weight':5}})
            self.lens_objects[standard_name]['z_Lens quality'].append({'value': 'spectroscopic', 'tracer': {'update status': 'in SILO', 'weight':5}})
            self.lens_objects[standard_name]['z_Source(s)'].append({'value': candidate['DETECTION_Z'], 'tracer': {'update status': 'in SILO', 'weight':5}})
            self.lens_objects[standard_name]['z_Source quality'].append({'value': 'spectroscopic', 'tracer': {'update status': 'in SILO', 'weight':5}})
            self.lens_objects[standard_name]['System Name'].append({'value': 'SDSS'+standard_name, 'tracer': {'update status': 'in SILO', 'weight':5}})
            if 'References' in self.lens_objects[standard_name]: self.lens_objects[standard_name]['References'].append('MNRAStmp(2021)303T')
            else: self.lens_objects[standard_name]['References'] = ['MNRAStmp(2021)303T']
            
                                                                                              
    def load_links(self):
        with open('LinKS_main.txt', 'r') as file:
            for line in file:
                line = line.split(',')
                if 'E-' in line[4]: line[4] = str(float(line[4]))
                if 'E-' in line[5]: line[5] = str(float(line[5]))
                standard_ra, standard_dec, standard_name = self.get_standard_name_and_coords({'RA [°]': line[4], 'Dec [°]': line[5]}, {'RA [°]':'RA [°]', 'Dec [°]':'Dec [°]'})
                if standard_name not in self.lens_objects: self.lens_objects[standard_name] = {'System Name':[], 'Discovery Date':[], 'RA (Hours part)':[], 'RA (Mins part)':[], 'RA (Secs part)':[], 'RA [°]': [], 'Dec (Degree part)': [], 'Dec (Arcmin part)': [], 'Dec (Arcsec part)': [], 'Dec [°]': [], 'Lens Grade': [], 'Number of images': [], 'Einstein_R ["]': [], 'z_Lens': [], 'z_Lens error':[], 'z_Lens quality':[], 'z_Source error':[], 'z_Source quality':[], 'z_Source(s)': [], 'Stellar velocity disp': [], 'Standard RA':[], 'Standard DEC':[], 'MLD_ID':[], 'Description':[], 'Lens type':[], 'Lens type MLD_ID':[], 'Discovery':[], 'Discovery_MLD_ID':[], 'MLD SDSS link':[], 'MLD ADSABS link':[], 'MLD NED link':[], 'MLD APOD link':[], 'References_MLD_ID':[], 'Has external link for SDSS':[], 'Has external link for ADSABS':[], 'Has external link for NED':[], 'Has external link for APOD':[]}
                self.lens_objects[standard_name]['Standard RA'].append({'value': standard_ra, 'tracer': {'update status': 'in LinKS', 'weight':0}})
                self.lens_objects[standard_name]['Standard DEC'].append({'value': standard_dec, 'tracer': {'update status': 'in LinKS', 'weight':0}})
                self.self.set_coord_details(standard_name, 5, 'LinKS', 'LinKS', 'in LinKS')
                self.lens_objects[standard_name]['System Name'].append({'value': line[0], 'tracer': {'update status': 'in LinKS', 'weight':0}})
                if 'References' in self.lens_objects[standard_name]: self.lens_objects[standard_name]['References'].append('MNRAS484(2019)3879P')
                else: self.lens_objects[standard_name]['References'] = ['MNRAS484(2019)3879P']
                    
        with open('LinKS_bonus.txt', 'r') as file:
            for line in file:
                line = line.split(',')
                standard_ra, standard_dec, standard_name = self.get_standard_name_and_coords({'RA [°]': line[1], 'Dec [°]': line[2]}, {'RA [°]':'RA [°]', 'Dec [°]':'Dec [°]'})
                if standard_name not in self.lens_objects: self.lens_objects[standard_name] = {'System Name':[], 'Discovery Date':[], 'RA (Hours part)':[], 'RA (Mins part)':[], 'RA (Secs part)':[], 'RA [°]': [], 'Dec (Degree part)': [], 'Dec (Arcmin part)': [], 'Dec (Arcsec part)': [], 'Dec [°]': [], 'Lens Grade': [], 'Number of images': [], 'Einstein_R ["]': [], 'z_Lens': [], 'z_Lens error':[], 'z_Lens quality':[], 'z_Source error':[], 'z_Source quality':[], 'Stellar velocity disp': [], 'Standard RA':[], 'Standard DEC':[], 'MLD_ID':[], 'Description':[], 'Lens type':[], 'Lens type MLD_ID':[], 'Discovery':[], 'Discovery_MLD_ID':[], 'MLD SDSS link':[], 'MLD ADSABS link':[], 'MLD NED link':[], 'MLD APOD link':[], 'References_MLD_ID':[], 'Has external link for SDSS':[], 'Has external link for ADSABS':[], 'Has external link for NED':[], 'Has external link for APOD':[]}
                self.lens_objects[standard_name]['Standard RA'].append({'value': standard_ra, 'tracer': {'update status': 'in LinKS', 'weight':0}})
                self.lens_objects[standard_name]['Standard DEC'].append({'value': standard_dec, 'tracer': {'update status': 'in LinKS', 'weight':0}})
                self.self.set_coord_details(standard_name, 5, 'LinKS', 'LinKS', 'in LinKS')
                self.lens_objects[standard_name]['System Name'].append({'value': line[0], 'tracer': {'update status': 'in LinKS', 'weight':0}})
                if 'References' in self.lens_objects[standard_name]: self.lens_objects[standard_name]['References'].append('MNRAS484(2019)3879P')
                else: self.lens_objects[standard_name]['References'] = ['MNRAS484(2019)3879P']
                    
    def convert_to_mld_reference_form(self, reference, url=None):
        """Prepare and process the bibcode to be converted from ads to MLD formatted bibcode"""
    
        if ('arXiv' in reference and '.' in reference) or 'FORJ0332' in reference:
            if url is None:
                reference = 'arXiv' + reference.split('arXiv')[-1]
                if ':' not in reference: reference = reference.replace('arXiv','arXiv:')
                url = join(self.base_links['ADS'], 'abs', reference)
            sleep(self.slow_down_seconds_after_requests)
            
            session = Session()
            response = session.get(url, headers = self.headers)
            check_response = self.check_passed(response)
            if check_response != 'Good': input('CHECK did not work%s'%check_response)
            r2 = response.text.split('<dt>Bibcode:</dt>\n                  <dd>\n                    <a href="/abs/')[1].split('/abstract">')[0]
            reference = r2
            session.close()
        if '(' not in reference: reference = self.convert_ads_to_mld_bibform(reference)
        return reference
                
        
    def convert_ads_to_mld_bibform(self, bibcode):
        """Convert the bibcode from ADS form to MLD form"""
        
        if 'arXiv' in bibcode and '.' in bibcode:
            year = bibcode[:4]
            bibparts = [bibpart for bibpart in bibcode.replace(year, '').split('.') if bibpart != '']
            mld_bibcode_form = bibparts[0] + '(' + year + ')' + bibparts[1]
            mld_bibcode_form = mld_bibcode_form[:-1] if mld_bibcode_form[-1] not in '1234567890' else mld_bibcode_form[:]
        elif 'arXiv' in bibcode: mld_bibcode_form = bibcode
        else:
            year = bibcode[:4]
            bibparts = [bibpart for bibpart in bibcode.replace(year, '').split('.') if bibpart != '']
            if len(bibparts) != 3: input('on hold to test %s %s'%(self.query, bibparts))
            mld_bibcode_form = bibparts[0] + bibparts[1] + '(' + year + ')' + bibparts[2]
            mld_bibcode_form = mld_bibcode_form[:-1] if mld_bibcode_form[-1] not in '1234567890' else mld_bibcode_form[:]
        return mld_bibcode_form 
        
    def check_referance_added_to_MLD_references(self, bibcode, title):
        """Check if either bibcode or title already in 'done' references"""
        
        mld_bibcode_form = self.convert_to_mld_reference_form(bibcode)
        if mld_bibcode_form in self.done_references_bib or title in self.done_references_title: return True
        else:
            self.done_references_bib.append(mld_bibcode_form)
            self.done_references_title.append(title)
            self.ads_to_mld_reference_interpreter[bibcode] = mld_bibcode_form
            self.update_reference.append(bibcode)
            return False
            
    def update_paper_overviews(self):
        """Set paper overviews"""
        
        self.set_papers_overview()
        for self.query in self.bibcodes:
            self.load_saved_overview()
            self.ads_scrapped_data[self.query]['Paper Overview'] = self.overview_list[self.query]
            self.save_overview()
            
    def set_papers_overview(self):
        """Load in the query information""" 
        
        converter = {'A1': 'first_authors', 'A2': 'secondary_authors', 'A3': 'tertiary_authors', 'A4': 'subsidiary_authors', 'AB': 'abstract', 'AD': 'author_address', 'AN': 'accession_number', 'AU': 'authors', 'C1': 'custom1', 'C2': 'custom2', 'C3': 'custom3', 'C4': 'custom4', 'C5': 'custom5', 'C6': 'custom6', 'C7': 'custom7', 'C8': 'custom8', 'CA': 'caption', 'CN': 'call_number', 'CY': 'place_published', 'DA': 'date', 'DB': 'name_of_database', 'DO': 'doi', 'DP': 'database_provider', 'EP': 'end_page', 'ER': 'end_of_reference', 'ET': 'edition', 'ID': 'id', 'IS': 'number', 'J2': 'alternate_title1', 'JA': 'alternate_title2', 'JF': 'alternate_title3', 'JO': 'journal_name', 'KW': 'keywords', 'L1': 'file_attachments1', 'L2': 'file_attachments2', 'L4': 'figure', 'LA': 'language', 'LB': 'label', 'M1': 'note', 'M3': 'type_of_work', 'N1': 'notes', 'N2': 'abstract', 'NV': 'number_of_Volumes', 'OP': 'original_publication', 'PB': 'publisher', 'PY': 'year', 'RI': 'reviewed_item', 'RN': 'research_notes', 'RP': 'reprint_edition', 'SE': 'version', 'SN': 'issn', 'SP': 'start_page', 'ST': 'short_title', 'T1': 'primary_title', 'T2': 'secondary_title', 'T3': 'tertiary_title', 'TA': 'translated_author', 'TI': 'Title', 'TT': 'translated_title', 'TY': 'type_of_reference', 'UK': 'unknown_tag', 'UR': 'url', 'VL': 'volume', 'Y1': 'publication_year', 'Y2': 'access_date'}
        
        #Scrapped IDs from MLD website and phrase from ADS reference abbreviations (http://cdsads.u-strasbg.fr/abs_doc/journals.html).
        journal_id_converter = {"0":"", "44":"Astrofizicheskie Issledovaniia Izvestiya Spetsial'noj Astrofizicheskoj Observatorii", "37":"arXiv e-prints", "11":"Astronomy and Astrophysics", "10":"The Astronomical Journal", "1":"The Astrophysical Journal", "8":"Astrophys. J. Lett.", "42":"The Astrophysical Journal Supplement Series", "45":"Bulletin d'information du telescope Canada-France-Hawaii", "2":"Monthly Notices of the Royal Astronomical Society", "39":"Nat", "38":"Nuc. Phy. B Proc. S.", "40":"Publications of the Astronomical Society of Japan", "43":"Revista Mexicana de Astronomia y Astrofisica", "41":"Science"}
         
        #Scrapped from MLD to ADS comparisons
        bibtype_converter = {"JOUR":"article", "Thesis/Dissertation":"phdthesis", "Preprint":"unpublished"}
        
        #Parse numbers missing in RIS version of RIS bib information from Leonidas list
        overview_number = {}
        with open(join(self.base_directory, 'resources', 'export-bibtexabs.bib'), 'r') as bibads_file:
            bibads_file = bibads_file.read()
            for line in bibads_file.split('\n'):
                if '@' in line: bibcode = line.split('{')[-1].split(',')[0]
                if 'number = ' in line: overview_number[bibcode] = line.split('{')[-1].split('}')[0]
            
        #Parse RIS formated bib information from Leonidas list
        with open(join(self.base_directory, 'resources', 'export-ris.txt'), 'r') as RIS_file:
            ris_file = RIS_file.read()
            delimiter = '  - '
            
            self.overview_list = {}
            for reference in ris_file.split('\n\n\n'):
                if reference == '': continue
                json_dictionary = {}
                for entry in reference.split('\n'):
                    if entry == 'ER  -': continue
                    key, value = entry.split(delimiter)
                    key = converter[key]
                    if key in ['authors', 'author_address', 'keywords']:
                        if key not in json_dictionary: json_dictionary[key] = []
                        json_dictionary[key].append(value)
                    else:
                        if key == 'url': bibcode = value.split('/')[-1]
                        if key in json_dictionary:
                            print('+++',key, json_dictionary.keys())
                            input('on hold')
                        json_dictionary[key] = value
                if bibcode in overview_number: json_dictionary['number'] = overview_number[bibcode]
                self.overview_list[bibcode] = json_dictionary
                

    def get_all_papers_referenced(self, references = [], name_versions = []):
        """Check all papers for references to this system"""
    
        for query in self.ads_scrapped_data.keys():
            total_text = ''
            if 'Paper text' in self.ads_scrapped_data[query]:
                total_text += self.ads_scrapped_data[query]['Paper text']
                if 'Table meta data' in self.ads_scrapped_data[query]: total_text += ''.join([self.ads_scrapped_data[query]['Table meta data'][table_set_key]['Response']for table_set_key in self.ads_scrapped_data[query]['Table meta data'].keys()])
            
            for name in name_versions:
                name2 = name.split('J')[-1]
                if name2 in total_text:
                    q2 = self.convert_to_mld_reference_form(query)
                    if q2 not in references:
                        references.append(self.convert_to_mld_reference_form(query))
                        break
        return references
            
    def parse_out_numbers(self, string):
        """Identifies numbers and decimals form a string and separates these in a vector. Used to strip either coordinates or error from column."""
        
        numbers = []
        temp = ''
        prev_s = ''
        if string == None: return []
        for s in string:
            if s in '0123456789.':
                if s == '.' and prev_s in '0123456789': temp+=s
                elif s in '0123456789': temp+=s
            else:
                if temp != '':
                    numbers.append(temp)
                    temp = ''
            prev_s = s*1
        if temp != '':
            numbers.append(temp)
            temp = '' #I know...not needed.
        return numbers
        
            
    def update_MLD_references(self):
        """Upload new references to Masterlens database"""
        #Trying session method failed, thus resorting to creating batch file for mysql update
        
        with open(join(self.base_directory, 'batch_update_mysql_references.txt'), 'w') as file:
            for self.query in self.update_reference:
                po = self.ads_scrapped_data[self.query]['Paper Overview']
                file.write("INSERT INTO reference ( siteID,identifier,abstract,author,title,journal,year,month,page,ads,keywords,editor,publisher,address,school,booktitle,series,bibtype,bibauthor,journalID,bibkey,public,modified")
                if 'number' in po: file.write(',number')
                if 'custom1' in po: file.write(',eprint')
                if "volume" in po: file.write(',volume')
                if "doi" in po: file.write(',doi')
                file.write(r' ) VALUES ( ')
                file.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,NOW()'%(1, self.ads_to_mld_reference_interpreter[self.query], po['abstract'], ', '.join(po['authors']), po['Title'], po['journal_name'], po['publication_year'].split('/')[0], po['start_page'], po['publication_year'].split('/')[1], po['start_page'], po['url'], ': '.join(po['keywords']), '', '', '', '', '', '', 'article', '', self.journal_id_converter_inverted[po['journalID']], 1))
                if 'number' in po: file.write(','+po['number'])
                if 'custom1' in po: file.write(','+po['custom1'].replace('eprint: ',''))
                if "volume" in po: file.write(','+po['volume'][0])
                if "doi" in po: file.write(','+po['doi'])
                file.write(' )\n')
            
    def update_MLD_lens_entries(self):
        """Update Masterlens database lens entries"""
        
        with open(join(self.base_directory, 'batch_update_mysql_lenses.txt'), 'w') as file:
            for system in self.lens_objects.keys():            
                add_system_dict = {"inputaction":"Save", "query_system_name": system}
                all_favoured_MLD = True
                none_favoured_in_MLD = True
                for key in self.lens_objects[system].keys():
                    if key not in self.masterlens_phrases_to_input_converter:
                        print('KEY>>>', key)
                        input('CHECK key should not be considered for upload')
                        continue
                    if len(self.lens_objects[system][key]) == 0: continue

                    weight = -9999
                    for data in self.lens_objects[system][key]:
                        if data['tracer']['weight'] >= weight:
                            if 'error' in data: value, method, error, weight, status = data['value'], data['method'], data['error'], data['tracer']['weight'], data['tracer']['update status']
                            else: value, method, error, weight, status = data['value'], '', '', data['tracer']['weight'], data['tracer']['update status']
                    if 'MLD' in status: none_favoured_in_MLD = False
                    elif 'MLD' not in status: all_favoured_MLD = False
                    methodid = '' if method in ['', 'MLD', 'NaN', ' NaN', None] else self.z_type_id[method] if 'z_lens' in key else self.z_type_id[method] if 'z_Source(s)' in key else self.er_quality_id[method] if 'Einstein_R ["]' in key else self.lens_type_id[method] if 'Lens type' in key else self.discovery_id[method] if 'Discovery' in key else ''

                    #THE MLD flag is a precaution for now that we can remove once a consensus on how we can update these correctly
                    if value == '': continue
                    add_system_dict[self.masterlens_phrases_to_input_converter[key]] = value
                    if (key + ' quality') in self.masterlens_phrases_to_input_converter: add_system_dict[self.masterlens_phrases_to_input_converter[key + ' quality']] = methodid
                    if (key + ' error') in self.masterlens_phrases_to_input_converter: add_system_dict[self.masterlens_phrases_to_input_converter[key + ' error']] = error

                add_system_dict['referencestoadd[]'] = '[' + ','.join([self.reference_id[reference] for reference in self.lens_objects[system]['References']]) + ']'
                add_system_dict['addreferences'] = 'addreferences'

                print('What POST looks like', add_system_dict)
                if all_favoured_MLD: 
                    print('Skipping since all entries in MLD for system:', system)
                    continue
                elif none_favoured_in_MLD:
                    file.write("INSERT INTO lens ( discovery_acronym,discovery_count,kind_acronym,kindID,filterID,system_name,lensgrade,multiplicity,morphology,reference_frame,equinox,description,alternate_name,z_lens,z_source,d_lens,d_source,vdisp,vdisp_err,time_delay0,time_delay1,mag_lens,mag_source,filter_lens,filter_source,theta_e,theta_e_err,theta_e_quality,theta_e_redshift,fluxes,ra_decimal,ra_hrs,ra_mins,ra_secs,ra_coord,ra_coord_err,dec_decimal,dec_degrees,dec_arcmin,dec_arcsec,dec_coord,dec_coord_err,number_images,reference_identifier,status,modified,created_by_member_name,modified_by_member_name,discovery_date,created,has_sdss,sdss_link,has_apod,apod_link,z_lens_err,z_lens_quality,z_source_err,z_source_quality,vett_status,released_status,hidden_status,vetted_by_member_name,released_as_of_version,released_by_member_name,hidden_by_member_name,vetted,released,hidden,repeats,graphic_status,coord_label,has_adsabs,adsabs_link,has_ned,ned_link,sdss_ObjID,sdss_specObjID,lens_name ) Values ( ")
                    file.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(self.discovery_id[add_system_dict['query_discoveryID_num1']], add_system_dict['query_discovery_count'], self.lens_type_id[add_system_dict['query_kindID']], add_system_dict['query_kindID'],0, add_system_dict['query_system_name'], add_system_dict['query_lensgrade'], '','','','J2000',  add_system_dict['query_description'], add_system_dict['query_alternate_name'],  add_system_dict['query_z_lens'], add_system_dict['query_z_source'], '','', add_system_dict['query_vdisp'],  add_system_dict['query_vdisp_err'], '','','','','','', add_system_dict['query_theta_e'],  add_system_dict['query_theta_e_err'], add_system_dict['query_theta_e_quality'], '','','', add_system_dict['query_ra_hrs'], add_system_dict['query_ra_mins'], add_system_dict['query_ra_secs'],  add_system_dict['query_ra_coord'], add_system_dict['ra_coord_err'], '', add_system_dict['query_dec_degrees'],  add_system_dict['query_dec_arcmin'], add_system_dict['query_dec_arcsec'], add_system_dict['query_dec_coord'], '','',referances[0],1,'NOW()', self.user_name, self.user_name,  add_system_dict['query_discovery_date'],'NOW()',None,None,None,None, add_system_dict['query_z_lens_err'],  add_system_dict['query_z_lens_quality'], add_system_dict['query_z_source_err'],  add_system_dict['query_z_source_quality'],0,1,0,'','','','','','','',0,0,'Manual',True,(r"https://ui.adsabs.harvard.edu/link_gateway/" + add_system_dict['referances'][0] + r"/PUB_HTML"),None,None,None,None, add_system_dict['query_system_name'].split('[')[0]))
                    file.write(' )\n')
                    print('SAVED SAVE NEW SYSTEM>>>>>', system)
                else: print('System in MLD but new information should be verified as to which entry to include...if any:', system)
