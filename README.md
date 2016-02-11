# Buckling-Check
#Buckling Check for DNV code with SNIF  
def buckling(a, b):
  #general input
  import math
  import numpy
  gamma_f = {'A':1.2, 'B':1.1}
  gamma_e = {'A':0.7, 'B':1.3}
  gamma_c = 1.0
  gamma_m = 1.15
  gamma_sc = 1.04
  gamma_epsilon = 2.0
  result ={}
  pc_root ={}
  lcc ={}
  dcc ={}
  m_sd ={}
  epsilon_sd={}
  s_sd ={}
  p_min = 0
  n = []
  g = []
  p_c = []
  f_y = a['alpha_u']*(a['SMYS']-a['f_ytemp'])
  f_u = a['alpha_u']*(a['SMTS']-a['f_utemp'])
  s_pt = f_y * 3.14159265 * (a['OD']- a['WT']) * a['WT']*0.001
  m_pt = f_y *(a['OD']- a['WT'])**2 * a['WT']*0.000001
  beta = (60 - a['OD']/a['WT'])/90
  alpha_c = (1- beta)+ beta * f_u/f_y
  p_el = (2 * a['E'] * (a['WT']/a['OD'])**3 )/(1-a['poisson_ratio']**2)
  p_p = f_y * a['alpha_fab'] * 2 * a['WT']/a['OD']
  if f_y < f_u/1.15:
    f_cb = f_y
  else:
    f_cb = f_u/1.15
  p_b = ((2* a['WT'])/(a['OD']-a['WT']))* f_cb * 2/(math.sqrt(3))
  print 'Node', 'Location','------DCC-----', '------LCC-----'
  for item in b:
    p_e = 0.01005525 * abs(b[item][7])
    #print p_e :check with excel
    epsilon_c = 100*0.78*((a['WT']/a['OD'])-0.01)*(1+5.75*((p_min-p_e)/p_b))*a['alpha_h']**(-1.5)*a['alpha_gw']
    #print epsilon_c, p_b:check with excel
    m = [1, -p_el, -1* (p_p**2+p_el*p_p*a['f_0']*(a['OD']/a['WT'])), p_el*p_p**2]
    f = numpy.roots(m)
    g = f.tolist()
    p_c = []
    for items in g:
      if items >= 0:
        p_c.append(items)
    m_sd ['A'] = gamma_f['A']* gamma_c * b[item][5]+ gamma_e['A']*(b[item][10]-b[item][5])
    m_sd ['B'] = gamma_f['B']* gamma_c * b[item][5]+ gamma_e['B']*(b[item][10]-b[item][5])
    epsilon_sd ['A'] = gamma_f['A']* gamma_c * abs(b[item][2]) + gamma_e['A'] * abs((abs (b[item][6])- abs (b[item][2])))
    epsilon_sd ['B'] = gamma_f['B']* gamma_c * abs(b[item][2]) + gamma_e['B'] * abs((abs (b[item][6])- abs (b[item][2])))
    s_sd ['A'] = gamma_f['A']* gamma_c * b[item][4]+gamma_e['A'] * ((b[item][9])- (b[item][4]))
    s_sd ['B'] = gamma_f['B']* gamma_c * b[item][4]+gamma_e['B']* ((b[item][9])- (b[item][4]))
    lcc ['A'] = (gamma_m * gamma_sc * (m_sd['A'])/(m_pt * alpha_c) + (gamma_m* gamma_sc * s_sd ['A'] /(alpha_c * s_pt))**2)**2+ (gamma_m* gamma_sc*(p_e-p_min)/(p_c[0]))**2
    lcc ['B'] = (gamma_m * gamma_sc * (m_sd['B'])/(m_pt * alpha_c) + (gamma_m* gamma_sc * s_sd ['B'] /(alpha_c * s_pt))**2)**2+ (gamma_m* gamma_sc*(p_e-p_min)/(p_c[0]))**2   
    dcc ['A'] = (epsilon_sd ['A']/(epsilon_c/gamma_epsilon))**0.8 + (p_e - p_min)/(p_c[0]/(gamma_m * gamma_sc))
    dcc ['B'] = (epsilon_sd ['B']/(epsilon_c/gamma_epsilon))**0.8 + (p_e - p_min)/(p_c[0]/(gamma_m * gamma_sc))
    result [item] = [round(dcc['A'], 3), round(dcc['B'], 3), round(lcc['A'],3) , round(lcc['B'],3)]
    print b[item][0], '  ' ,b[item][1], result[item]
  return
 
 
def pipe_properties(filename):
# working well, gives all buckling check as a dict which was saved on csv file
  import numpy
  with open (filename, 'r') as myfile:
    data = myfile.read()
    my_data= data.split('\n')
    parameters ={}
    for item in my_data:
      raw = item.split(',')
      parameters [raw[0]] = float(raw[1])
  #for item in parameters.keys():
   #parameters[item]=numpy.array(parameters[item]) 
  #print parameters, type(parameters)
  return parameters  

def data(filename):
  import numpy
  m = {}
  n = {}
  f = {}
  h = {}
  final ={}
  final2 ={}
  with open (filename, 'r') as myfile:
    data0 = myfile.read()
    data = data0.replace ("'" , "")
    line_data= data.split('\n')
    for i in range(0, (len(line_data)/2)):
      m [i]= 'static: '+ line_data[i] + ' dynamic: ' + line_data [i + len(line_data)/2]
  for item in m.keys():
    #print m[item]
    n[item] = m[item].replace ('[', ',')
    f[item] = n[item].replace (']', ',')
    h[item] = f[item].replace (',', ' ')
  for item in h.keys():
    final[item] = []
  for item in h.keys():
    for i in range (1, 15):
      if i != 7 and i != 8 and i != 9: 
        final[item].append(h[item].split()[i])
  #for item in final.keys():
    #print final[item]
  for thing in final.keys():
    for item in final[thing][2:]:
      flt = float(item)
      g = final[thing].index (item)
      final[thing][g] = flt 
  for item in final.keys():
    final[item][0] = int(final[item][0])
    #print item, final[item]
  return final
       
import operator
import random
import re
import os
import sys

def main():    
  a = pipe_properties('buckling paramters.csv')
  b = data ('dy_st.OUT')
  buckling (a, b)
    
if __name__ == '__main__':
  main()
