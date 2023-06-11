import numpy as np
from dataclasses import dataclass          #required to declare a type

from astropy.io import fits as fits
from astropy.time import Time

#hdu.header['USERNAME'] = "Mao"
#hdu.header["INSTRU"]="GMG"
#hdu.header["OBSDATE"]="2022-04-12T12-40-37"
#hdu.header['FILTER'] = "sdssr"
#hdu.header['TARGET'] = "GRB220412"
#hdu.header['STACK'] = 1
#hdu.header['GAIN'] = 0.33
#hdu.header["EXP"]=900

@dataclass
class OBS:  
	instru : str
	obserdate : str
	mjd : float
	exp: str
	username: str
	filters: str
	target: str
	stack: str
	gain: float
	filter_ref:str
	magsys: str
	catalogcalib: str
	filtercalib: str
	catalogref: str

	def __init__(self):
		self.instru = "None"
		self.obserdate = "None"
		self.mjd= 0.0
		self.exp = "None"
		self.username = "None"	
		self.filters = "None"	
		self.target = "None"	
		self.stack = "None"	
		self.gain = 1.0
		self.magsys="AB"
		self.catalogcalib="ps1"
		self.filtercalib="rmag"
		self.catalogref="PanSTARRS/DR1/r"
		self.filter_ref="rmag"

def test_instrument(instrument,valid):
	GRANDMA={'SBO','SRO','TCA','TCH','TRE','FZU-Auger','FZU-CTA-N','Xinglong-2.16m','TNT','SNOVA','ALi-50','OPD-1.6m','OPD-60','SOAR','Abastumani-T70','Abastumani-T48','C2PU','GMG','UBAI-T60N','UBAI-T60S','Makes-T60','MOSS','OWL','TRAPPIST','ShAO-T60','ShAO-2m','CAHA','VIRT','OST','PicDuMidi-T1M','CFHT','GTC','TELESCO','KAO','NOWT','Xinglong-2.16m','OST-CDK','Lisnyky','ShAOT60','UBAI-ST60','CFHT-Megacam',"ASTEP"}

	KNC={'C11FREE','HAN','K26','SUTO1','SUTO2','T-BRO','T-CAT','TJML','TJMS-PLANETE-SCIENCES','T-AGU','T19','HAO','T24','T21'}
	if (instrument in KNC) or (instrument in GRANDMA): 
		return valid
	else:
		print("wrong name instrument")
		return 0
		
def 	test_username(username,valid):
	GRANDMA={'Noysena','Simon','Klotz','DeUgarte','Karpov','Zhu','Letian','Masek','Blazek','Zhou','Du','Zeng','Song','Corradi',
	'Navarete','Kochiashivili','Beradze','Rivet','Mao','Burkhonov','Benkhaldoun','Gurbanov','Hesenov','Vidadi','Agayeva','Ismailovn',
	'Kann','Orange','Gokuldass','Hainich','Takey','Fouad','Elhosseiny','Hellmich','Rachith','Irureta','Wangletian','Rinner','Guo','Songxuan','Emin-emrah','Antier',"Pormente","Abastumani-T70"}		
	KNC={'Boust','Boutigny','Vallieres','Broens','Cailleau','Cejudo','Eggenstein','Bertrand','Bayard','Menard','Freeberg','Galdies',
	'Granier','Jaquiery','Leonini','Leroy','Popowicz','Rousselot','Serrau','Kneip','Richmond','Oksanen','Marchais','Aguerre','Kaeouech','Rinner-Benkhaldoun'}
	if (username in KNC) or (username in GRANDMA):
		return valid
	else:
		print("wrong name username")
		return 0

def 	test_filter(filters,valid):
	#filters='sdssr'
	if filters=="SDSS-z'":
		filters='sdssz'
	if filters=="SDSS-r'":
		filters='sdssr'
	if filters=="SDSS-r'":
		filters='sdssr'
	if filters=="SDSSr":
		filters='sdssg'
	if filters=="SDSS-i'":
		filters='sdssi'
	if filters=="R/Bessel":
		filters="R-Bessel"
	if filters=="Red":
		filters="R"
	if filters=="SDSS-g'":
		filters="sdssg"	
	if filters=="g.MP9402":
		filters='sdssg'	
	if filters=="r.MP9602":
		filters='sdssr'	
	if filters=="i.MP9703":
		filters="sdssi"
	if filters=="z.MP9901":
		filters="sdssz"
	if filters=="sdrssr":
		filters="sdssr"		

	if filters in {'I','Ic','C','Clear','clear','L','lumen','LUM','lum','G','TG','g','sG','sdssi','sloanG','Green','R','Rc','r','sR','sloanR','SR','g','sG','i','sloanI','v','V','TB','TR','R-Bessel','bessellr','SDSSrp+','sdssr','sdssz','sdssg',"sdssi","SDSS-g","RGaia"}:
		return valid,filters
	else:
		print("wrong name filter")
		return 0

def get_gain(instru):
	if instru=="KAO":
		return 2.14
	if instru=="Lisnyky":
		return 0.82
	if instru=="MOSS":
		return 100.
	if instru=="ShAOT60":
		return 1.250
	if instru=="UBAI-ST60":
		return 3.2
	if instru=="UBAI-AZT-22":
		return 1.45
	if instru=="SNOVA":
		return 2.0
	if instru=="Abastumani-T70":
		return 1.35
	if instru=="TCH":
		return 0.91
	if instru=="HAO":
		return 200.
	if instru=="CFHT-Megacam":
		#return 1.5611
		#return 1.6418 
		return 1.6640 
	if instru=="GTC":
		#return 1.5611
		#return 1.6418 
		return 0.95
def assess_refcat(obs):
	if obs.filters=="sdssz":
		if obs.instru=="KAO" or obs.instru=="CFHT-Megacam":
			obs.catalogcalib="ps1"
			obs.filtercalib="zmag"
			obs.filter_ref="zmag"
			obs.catalogref="PanSTARRS/DR1/z"
	if obs.filters=="sdssr":			
		if obs.instru=="KAO" or obs.instru=="CFHT-Megacam":
			obs.catalogcalib="ps1"
			obs.filtercalib="rmag"
			obs.filter_ref="rmag"
			obs.catalogref="PanSTARRS/DR1/r"
		if obs.instru=="GTC":
			obs.catalogcalib="ps1"
			obs.filtercalib="rmag"
			obs.filter_ref="rmag"
			obs.catalogref="PanSTARRS/DR1/r"
	if obs.filters=="Clear":			
		if obs.instru=="SNOVA":
			obs.catalogcalib="ps1"
			obs.filtercalib="rmag"
			obs.filter_ref="rmag"
			obs.catalogref="PanSTARRS/DR1/r"
	if obs.filters=="R/Bessel":			
		if obs.instru=="UBAI-ST60":
			obs.catalogcalib="ps1"
			obs.filtercalib="rmag"
			obs.filter_ref="rmag"
			obs.catalogref="PanSTARRS/DR1/r"
	if obs.filters=="R":			
		if obs.instru=="Lisnyky":
			obs.catalogcalib="ps1"
			obs.filtercalib="rmag"
			obs.filter_ref="rmag"
			obs.catalogref="PanSTARRS/DR1/r"		
		if obs.instru=="ShAOT60":
			obs.catalogcalib="ps1"
			obs.filtercalib="rmag"
			obs.filter_ref="rmag"
			obs.catalogref="PanSTARRS/DR1/r"
	if obs.filters=="I":			
		if obs.instru=="Abastumani-T70":
			obs.catalogcalib="ps1"
			obs.filtercalib="imag"
			obs.filter_ref="imag"
			obs.catalogref="PanSTARRS/DR1/i"
	if obs.filters=="sdssg":			
		if obs.instru=="KAO" or obs.instru=="CFHT-Megacam":
			obs.catalogcalib="ps1"
			obs.filtercalib="gmag"
			obs.filter_ref="gmag"
			obs.catalogref="PanSTARRS/DR1/g"
	if obs.filters=="sdssi":			
		if obs.instru=="KAO" or obs.instru=="CFHT-Megacam":
			obs.catalogcalib="ps1"
			obs.filtercalib="imag"
			obs.filter_ref="imag"
			obs.catalogref="PanSTARRS/DR1/i"
	if obs.filters=="I":			
		if obs.instru=="Abastumani":
			obs.catalogcalib="ps1"
			obs.filtercalib="imag"
			obs.filter_ref="imag"
			obs.catalogref="PanSTARRS/DR1/i"
	if obs.filters=="RGaia":			
		if obs.instru=="ASTEP":
			obs.catalogcalib="gaiaedr3"
			obs.filtercalib="RPmag"
			obs.filter_ref="BPmag"
			obs.catalogref='flux-Rp/I/350/gaiaedr3'
	if obs.filters=="Green":			
		if obs.instru=="ShAOT60":
			obs.catalogcalib="ps1"
			obs.filtercalib="gmag"
			obs.filter_ref="gmag"
			obs.catalogref="PanSTARRS/DR1/g"
	return obs

def Cformatandcollectinfo(filename):

	data_fits=""
	
	data_fits = [filename]
	
	#print(data_fits)
	obs: OBS = OBS()
	for filefits in data_fits:
		print(filefits.split('_'))
		valid=1
		valid2=1
		hdu = fits.open(filefits)
		hdr=hdu[0].header		
		#USERNAME: follow strictly the list attached
		#INSTRU: follow strictly the list attached
		#OBSDATE: start of the observation (the date when the first image has been taken) in UTC
		#FILTER: follow strictly the list
		#TARGET: XX
		#STACK: (0,1,...N)
		#GAIN: gain
		#EXP: exposure time in the format of 10x100s
		addusername=""
		addinstru=""
		addobsdate=""
		addfilter=""
		addtarget=""
		addstack=""
		addgain=""
		addexp=""
		instru=""
		obserdate=""
		try:
						obs.username=hdr["USERNAME"]
						print("username",obs.username)
		except KeyError:
			valid=0
			addusername="#USERNAME keyword,"
			try:
				obs.username = filefits.split('_')[2].capitalize()
				valid2=test_username(obs.username,valid2)
				#print(valid2)
			except:
				print("wrong filename")		
			
		try:
						obs.instru=hdr["INSTRU"]
						if obs.instru=="C2PU/Omicron":
							obs.instru="C2PU-Omicron"
						print("instru",obs.instru)
		except KeyError:
			valid=0
			addinstru="#INSTRU keyword,"
			try:
				obs.instru = filefits.split('_')[3]
				valid2=test_instrument(obs.instru,valid2)
			except:
				print("wrong filename")		
		try:
						date=hdr["OBSDATE"]
						obs.obserdate=date.split('T')[0]+'T'+date.split('T')[1].replace('-',':')
						t = Time(obs.obserdate, format='isot', scale='utc')
						obs.mjd=t.mjd
						print("obsdate",obs.obserdate)
		except KeyError:
			valid=0
			addobsdate="#OBSDATE keyword (start of the obs)"
			try:
				date = filefits.split('_')[4]
				obs.obserdate=date.split('T')[0]+'T'+date.split('T')[1].replace('-',':')
				t = Time(obs.obserdate, format='isot', scale='utc')
				obs.mjd=t.mjd
			except:
				print("wrong filename")				
		try:
						obs.filters=str(hdr["FILTER"])
						#print(obs.filters)
						valid2,obs.filters=test_filter(obs.filters,valid2)
						obs=assess_refcat(obs)
						print("filters",obs.filters)
		except KeyError:
			valid=0
			addfilter="#FILTER keyword,"
			try:
				obs.filters = filefits.split('_')[5]
				valid2,obs.filters=test_filter(obs.filters,valid2)	
				obs=assess_refcat(obs)
			except:
				print("wrong filename")		
		try:
						obs.target=hdr["TARGET"]
						print("target",obs.target)
		except KeyError:
			valid=0
			addtarget="#TARGET keyword, "	

		try:
						obs.stack=hdr["STACK"]
						print("stack",obs.stack)
		except KeyError:
			valid=0
			addstack="#STACK keyword, "		

		try:
						obs.gain=float(hdr["GAIN"])
						print("gain",obs.gain)
		except KeyError:
			valid=0
			addstack="#GAIN keyword, "	
			obs.gain=get_gain(obs.instru)
		
		try:
						exp=hdr["EXPO"]
						print("exp",obs.exp)
		except KeyError:
			valid=0
			addex="#EXP keyword, "	
			try:
				obs.exp = filefits.split('_')[7].split(".")[0]
				#print(filefits.split('_')[7])
			except:
				print("wrong filename")		
		#print(valid)							
		if valid ==0:
			print("Format problem in "+str(filefits), "inexistant keywords: "+addusername+addinstru+addobsdate+addfilter+addtarget+addstack+addexp)
	return obs

def treat_GRBdata(filename):
		obs=Cformatandcollectinfo(filename)
		return obs

		
