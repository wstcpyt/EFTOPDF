#read NORMAL INCIDENCE E-field profiles and calculate absorption and excitong generation, PDF per wavelength
import numpy as np
from numpy import *
import pylab


#constants
h=6.6260755E-34         #planks
c=299792458             #speed of light
epsilon=8.85418782E-12  #permittivity in vacuum
e=1.60217733E-19        #charge of electron
pangle=52.427           #prism angle (52.427 + (asin(sin(col("1")*3.14159/180)/2.1868)*180/3.14159))
aerror=0.3              #difference in the angle between data and baseline
bcurrent=1E-9           #baseline current (find from the plot)
criticalangle=-30.0     #critical angle
ecutoffpos=-0.01        #electron cutoff position for 10percent value
hcutoffpos=-0.01        #hole cutoff  position or 10percent value

#photoenergy
def photonenergy(wl):
    return(h*c/(wl*1E-9))

#calculate absorption rate
def absorption(efield,n,k,energy):
     return((2*np.pi*epsilon*k*n*efield)*1E-6/h)

#optical constants
Mo_n = np.genfromtxt('oc/Mo-n.txt', delimiter='\t',invalid_raise=False)        #Au n
Mo_k = np.genfromtxt('oc/Mo-k.txt', delimiter='\t',invalid_raise=False)        #Au k
CdS_n = np.genfromtxt('oc/CdS-n.txt', delimiter='\t',invalid_raise=False)    #v2o5 n
CdS_k = np.genfromtxt('oc/CdS-k.txt', delimiter='\t',invalid_raise=False)    #v205 k
azo_n = np.genfromtxt('oc/ZnO-n.txt', delimiter='\t',invalid_raise=False)      #ITO n
azo_k = np.genfromtxt('oc/ZnO-k.txt', delimiter='\t',invalid_raise=False)      #ITO k
zno_n = np.genfromtxt('oc/ZnO-n.txt', delimiter='\t',invalid_raise=False)    #TiO2 n
zno_k = np.genfromtxt('oc/ZnO-k.txt', delimiter='\t',invalid_raise=False)    #TiO2 k
CIGS_n_ipl = np.genfromtxt('oc/CIGS-Ga-p3-n.txt', delimiter='\t',invalid_raise=False)    #CIGS-n-inplane
CIGS_k_ipl = np.genfromtxt('oc/CIGS-Ga-p3-k.txt', delimiter='\t',invalid_raise=False)    #CIGS-k-inplane
CIGS_n_opl = np.genfromtxt('oc/CIGS-Ga-p3-n.txt', delimiter='\t',invalid_raise=False)    #CIGS-n-outplane
CIGS_k_opl = np.genfromtxt('oc/CIGS-Ga-p3-k.txt', delimiter='\t',invalid_raise=False)    #CIGS-k-outplane


#get n and k
def nandk(material,wl):
    if material=='Mo':
        return (Mo_n[list(Mo_n[:,0]).index(wl),1],Mo_k[list(Mo_k[:,0]).index(wl),1]) #convert array to list,search wl,get index,
    if material=='CdS':
        return (CdS_n[list(CdS_n[:,0]).index(wl),1],CdS_k[list(CdS_k[:,0]).index(wl),1])
    if material=='CIGS': #for TE
        return (CIGS_n_ipl[list(CIGS_n_ipl[:,0]).index(wl),1],CIGS_k_ipl[list(CIGS_k_ipl[:,0]).index(wl),1])
    if material=='CIGS_opl':
        return (CIGS_n_opl[list(CIGS_n_opl[:,0]).index(wl),1],CIGS_k_opl[list(CIGS_k_opl[:,0]).index(wl),1])
    if material=='ZnO':
        return (zno_n[list(zno_n[:,0]).index(wl),1],zno_k[list(zno_k[:,0]).index(wl),1])
    if material=='AZO':
        return (azo_n[list(azo_n[:,0]).index(wl),1],azo_k[list(azo_k[:,0]).index(wl),1])


#get pos,value arrays out from the internal data structures
def get_internal(odata):
    pos   = zeros(odata.size)
    val   = zeros(odata.size)
    for i in range(0,odata.size,1):
        pos[i],val[i]=odata[i][0],odata[i][1]
    return pos,val



#define CIGS structure
dt = np.dtype([('material', np.str_, 16), ('position', np.float64, (2,))])
device=np.array([('AZO',(-0.1,-0.275)),              #AZO
                 ('ZnO',(-0.275,-0.365)),                #ZnO
                 ('CdS',(-0.365,-0.415)),            #exact
                 ('CIGS',(-0.415,-2.715)),       #CIGS-0p3
                 ('Mo',(-2.715,-3.045))],dtype=dt)   #Mo
res=0.001 #resolution


#find material given position
def material(pos,device):
    for i in device:
        if abs(pos) >= abs(i['position'][0])and abs(pos) <= abs(i['position'][1]):
            return (i[0])

#savefile
def savedata(poutdata):
    out=np.zeros((poutdata.size,7)) #wl,mode,abs,peak_position,mean_position,10nm_tio2,10nm_v2o5
    j=0;
    for row in poutdata:
        out[j]=(row[0],
                row['mode_no'],
                row['absorption'],
                row['peak_position'],
                row['mean_position'],
                row['10nm_from_tio2_exciton_percentage'],
                row['10nm_from_v2o5_exciton_percentage'])
        j=j+1
    return out

#get normalized position in the active layer##########################
def normalizepos(pos):
    js=list(device['material']).index('CdS')
    je=list(device['material']).index('Mo')
    return((pos-device[js][1][1])/(device[je][1][0]-device[js][1][1]))
#####################################################################

#get normalized position in the active layer##########################
def normalizefullpos(pos):
    js=list(device['material']).index('AZO')
    je=list(device['material']).index('Mo')
    return((pos-device[js][1][0])/(device[je][1][1]-device[js][1][0]))
#####################################################################


#Parameter file giving the wavelengths, te mode # and tm mode #
parafilename='wavelength_normal.txt'
parafile = np.genfromtxt(parafilename, delimiter='\t',invalid_raise=False,dtype=np.uint16)

#define output data structure
odt1 = np.dtype([('position', np.float64,1), ('abs_rate', np.float64,1)])                   #for reference not used
gensize = np.absolute(np.round_((device[4]['position'][1]-device[0]['position'][0])/res))+1 #size of the generation array
pdfsize = np.absolute(np.round_((device[4]['position'][0]-device[2]['position'][1])/res))   #size of the propability distribution array
outputdt = np.dtype([('wavelength', np.uint16,1),
                     ('mode',       np.str_,1),
                     ('mode_no',    np.uint8,1),
                     ('beta_fdtd',  np.float64,1),
                     ('absorption',  np.float64,1),
                     ('generation_rate',np.dtype([('position', np.float64,1), ('abs_rate', np.float64,1)]),(gensize,)),
                     ('active_gen_rate',np.dtype([('position', np.float64,1), ('abs_rate', np.float64,1)]),(pdfsize,)),
                     ('pdf',            np.dtype([('position', np.float64,1), ('pdf', np.float64,1)]),(gensize,)), #pdfsize
                     ('cdf',            np.dtype([('position', np.float64,1), ('cdf', np.float64,1)]),(gensize,)), #pdfsize
                     ('peak_position', np.float64,1),
                     ('mean_position', np.float64,1),
                     ('10nm_from_CdS_exciton_percentage', np.float64,1),
                     ('10nm_from_Mo_exciton_percentage', np.float64,1)
                     ])
outdata = np.zeros(parafile.size,dtype=outputdt)   #output file (#ofmodes)

#make plot generation_rate
fig1= pylab.figure()
slegend1=np.array([])
p=fig1.add_subplot(111)
#p.set_title("NORMAL INCIDENCE")
p.set_ylabel('NORMAL INCIDENCE Generation PDF (um-1)')#('Generation Rate (m-3Ss-1)')
p.set_ylabel('PDF(1/um)')#('Generation Rate (m-3Ss-1)')
p.set_xlabel('$x$/$d$')
#pylab.grid(True)
#linetype=array([('b-','bo'),('g-','go'),('r-','ro'),('c-','co'),
#                ('m-','mo'),('y-','yo'),('k-','ko'),
#                ('bx','b*'),('gx','g*'),('rx','r*'),('cx','c*'),
#                ('b.','b--'),('g.','g--'),('r.','r--'),('c.','c--')])
#linetype=array([('b-','b--'),('g-','g--'),('r-','r--'),('c-','c--'),
#                ('m-','m--'),('y-','y--'),('k-','k--'),('bx','b--'),
#                ('gx','gt'),('rx','rt'),('cx','ct'),
#                ('b<','bt'),('g<','go'),('r<','ro'),('c<','co')])
linetype=array([('b-','b--','b-.','b:'),('g-','g--','g-.','g:'),('r-','r--','r-.','r:'),('c-','c--','c-.','c:'),
                ('m-','m--','m-.','m:'),('y-','y--','y-.','y:'),('k-','k--','k-.','k:'),('b-.','bx--','bs.','b>'),('g-.','gx','gs.','g>'),
                ('r-.','rx','rs','r>'),('b--','b-','b-.','b:'),('g--','g-','g-.','g:'),('r--','r-','r-.','r:'),('c--','c-','c-.','c:'),
                ('m--','m-','m-.','m:'),('y--','y-','y-.','y:'),('k--','k-','k-.','k:')])
pylab.rcParams.update({'font.size': 12})
pylab.rcParams.update({'font.name': 'Arial'})



###############################################################
#calculate abs,generatio rate, pdf for each wavelength and mode
print('reading parameter file:%s'%parafilename)
j=ln=0
for para in parafile:
    wl =  para #wavelength

    outdata[j]['wavelength']= wl
    outdata[j]['mode']= 'N'
    outdata[j]['mode_no']= 0
    print('Reading.. CIGS-EF-Ga-p3-%dnm-Ex'%(wl))
    df = np.genfromtxt('data/CIGS-EF-Ga-p3-%dnm-Ex.txt'%(wl),delimiter='\t',skip_header=955,invalid_raise=False) #DELETED THE GLASS REGION
    #exciton generation rate-total###################################
    i=0
    for field in df:
        #print(field[0])
        #print(material(field[0],device))
        n,k = nandk(material(field[0],device),wl)
        outdata[j]['generation_rate'][i]= array([field[0],absorption(field[1],n,k,photonenergy(wl))])
        i=i+1
    print('\tAnalysis wl:%d,Normal'%(wl))
    #pos,val=get_internal(outdata[j]['generation_rate'])
    #slegend1=np.append(slegend1,'%d,NORMAL'%(wl),axis=None)
    #p.plot(normalizefullpos(pos),val,(linetype[ln,0]),linewidth=2.0) #pick line color from array
    #tempos=zeros(pos.size)
    #tempos.resize(pos.size,2)
    #tempos[:,0]=normalizefullpos(pos)
    #tempos[:,1]=val
    #np.savetxt('%d_NORMAL_GENERATION_AU.dat'%(wl),tempos,delimiter='\t')

    #absorption-active layer##########################################
    tsum=tactive=0
    for field in df:
        n,k = nandk(material(field[0],device),wl)
        tsum=tsum + absorption(field[1],n,k,photonenergy(wl))
        if material(field[0],device) == 'CIGS':
            tactive=tactive+absorption(field[1],n,k,photonenergy(wl))
    outdata[j]['absorption'] = tactive/tsum
    print('\tAbsorption wl:%dnm:NORMAL:%f...'%(wl,outdata[j]['absorption']))

    #exciton generation rate-active layer##############################
    i=0
    for field in df:
        if material(field[0],device) == 'CIGS':
            n,k = nandk(material(field[0],device),wl)
            outdata[j]['active_gen_rate'][i]= array([field[0],absorption(field[1],n,k,photonenergy(wl))])
            i=i+1

    pos,val=get_internal(outdata[j]['active_gen_rate'])
    slegend1=np.append(slegend1,'%d,NORMAL'%(wl),axis=None)
    #p.plot(normalizepos(pos),val,(linetype[ln,0]))


    #statistical analysis of the generation profile####################
    pos1,val=get_internal(outdata[j]['generation_rate'])    #'active_gen_rate'
    pos=normalizefullpos(pos1)  #normalizepos
    tgeneration=np.trapz(val,x=pos) #integration (the x axis is in um's)
    outdata[j]['cdf'][0]=array([pos[0],0]) #initialize cumalative distribution
    for i in range(0,pos.size,1):
        outdata[j]['pdf'][i]= array([pos[i],(val[i]/np.absolute(tgeneration))]) #pdf
        if i>0:
            a=outdata[j]['pdf'][i-1][1]
            b=outdata[j]['pdf'][i][1]
            last=outdata[j]['cdf'][i-1][1]
            outdata[j]['cdf'][i]= array([pos[i],(((a+b)*abs(pos[i-1]-pos[i])*0.5)+last)]) #cdf

    #check if the PDF is 1
    pos,val=get_internal(outdata[j]['pdf'])
    tgeneration=np.trapz(val,x=pos)
    print('\tCDF of %dnm:%f'%(wl,tgeneration))

    pos,val=get_internal(outdata[j]['pdf'])
    mean=np.sum(array(val[:]*pos[:]))/np.sum(val)  #calculating mean gen position
    outdata[j]['mean_position']=mean #normalizepos(mean) #mean
    outdata[j]['peak_position']=pos[list(val).index(val.max())] #normalizepos(pos[list(val).index(val.max())]) #peak
    print('\tmean:%f\n\tpeak:%f..'%(outdata[j]['mean_position'],outdata[j]['peak_position']))
    '''
    #exciton generation fraction at x nm from interface
    fromtio21=(device[list(device['material']).index('CdS')]['position'][0])-ecutoffpos #from TiO2
    fromtio2=normalizepos(fromtio21)
    pos,val=get_internal(outdata[j]['cdf'])
    outdata[j]['10nm_from_CdS_exciton_percentage']=val[list(pos).index(fromtio2)]
    print('\texcitons within:%d from CdS:%f..'%(ecutoffpos,outdata[j]['10nm_from_CdS_exciton_percentage']))

    fromv2o51=(device[list(device['material']).index('Mo')]['position'][1])+hcutoffpos #from V2O5
    fromv2o5=normalizepos(fromv2o51)
    pos,val=get_internal(outdata[j]['cdf'])
    outdata[j]['10nm_from_Mo_exciton_percentage']=val[val.size-1]-(val[list(pos).index(fromv2o5)])
    print('\texcitons within:%d from Mo:%f..'%(hcutoffpos,outdata[j]['10nm_from_Mo_exciton_percentage']))
       '''
    #pos,val=get_internal(outdata[j]['pdf'])
    #slegend1=np.append(slegend1,'%d,NORMAL'%(wl),axis=None)
    #p.plot(pos,val,linetype[ln,0],linewidth=2.0) #pick line color from array

    np.savetxt('%d_NORMAL_PDF.dat'%(wl),outdata[j]['pdf'],delimiter='\t')


    j = j+1 #next wl:mode

    ln=ln+1 #forplotcolor


#format legend for 'generation_rate'
leg1=p.legend(slegend1, 'upper right', shadow=True)
for t in leg1.get_texts():
    t.set_fontsize('small')


#save data
#np.savetxt('NORMAL-mode_abs_peak_mean_fromtio2_fromv2o5-SHORT-AU.txt',savedata(outdata),delimiter=',') #save wl,mode_no,abs
pylab.show()


