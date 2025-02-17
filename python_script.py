import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os


'''
  Fix Earth's declination
'''
beta = 23.5 * np.pi/180 # declication



def solve_cos_trig(A):
  '''
    Solve equation cos(x-phi)>=A
    and return the Lebesgue measure (length) of the set os solution limited to [0, 2*pi]
    phi is in range -pi/2 to pi/2

    The solution turns out to be very simple
  '''
  if A >= 1:
    return 0
  if A <=-1:
    return 2*np.pi
  else:
    return 2*np.arccos(A)


def solve_sin_trig(A):
  '''
    Solve equation sin(x)>=A
    and return the Lebesgue measure (length) of the set os solution limited to [0, 2*pi]

    The solution turns out to be very simple
  '''

  if A >= 1:
    return 0
  if A <=-1:
    return 2*np.pi
  if A > 0:
    return np.pi - 2 * np.arcsin(A)


def solve_trig(A,B,C):
  '''
    Solves equation A * cos(x) + B * sin(x) >= C
    and returns the Lebesgue measure (length) of the set of solutions limited to [0, 2*pi]
  '''
  if A == 0 and B == 0:
    if C <= 0:
      return 2*np.pi
    if C > 0:
      return 0

  if A == 0 and B != 0:
    return solve_sin_trig(C/B)

  if A !=0 and B == 0:
    return solve_cos_trig(C/A)

  if A!=0 and B!=0:
    return solve_cos_trig(C/(A**2+B**2)**0.5)

  return None



def plot_daylight(gamma, df=None, num_years=1, title=None, vis=False):
  '''
    Plots the estimated daylight for a full year

    If df - pandas data frame is passed, it will also plot the estimated
    daylight times with the ones in the df. In this case df must containg columns
    "day" and "daylight"
  '''

  res = []

  nine_days = (9/365) * 2*np.pi # t=0 corresponds to the 22nd of Decemver. Shift required
  sun_shift =   - 0.016675632060247445 # how much sun is moved relative to the center of the orbit
  for t in np.arange(nine_days, (num_years * 2*np.pi + nine_days), 2*np.pi/365):
    A = np.cos(gamma)*np.cos(beta)*(np.cos(t) + sun_shift)
    B = np.cos(gamma)*np.sin(t)
    C = np.sin(beta)*np.sin(gamma)*(np.cos(t) + sun_shift)
    res.append(24 * solve_trig(A,B,C)/(2*np.pi))

  ok = False # will be true if proper df is passed
  rmse = None
  max_diff = None

  if df is not None:
    ok = 'day' in df.columns and 'daylight' in df.columns

  if ok:
      plt.plot(df['daylight'], color = 'red', label = 'daylight from file')
      rmse = np.sqrt(np.mean((df['daylight'] - res)**2))
      max_diff = np.max(np.abs(df['daylight'] - res))

  if vis:
    plt.plot(res, color = 'blue', label='estimated daylight')
    plt.ylim(0, 24.5)
    plt.xlabel('day')
    plt.ylabel('daylight [h]')
    plt.legend()
    plt.title(title)
    plt.show()

  return {'result': res, 'rmse': rmse, 'max_diff': max_diff }

'''
  Basic checks
   - check if daylight is 12h constant on equator
   - check if daylight is 0h or 24h on north and south poles
'''
a = plot_daylight(0, df=None, num_years=1, title='Daylight on equator')
b = plot_daylight(np.pi/2, df=None, num_years=1, title='Daylight on north pole')
c = plot_daylight(-np.pi/2, df=None, num_years=1, title='Daylight on south pole')
d = plot_daylight(75*np.pi/180, df=None, num_years=1, title='Above the arctic circle')
e = plot_daylight(10*np.pi/180, df=None, num_years=1, title='Close to equator')


fig, ax = plt.subplots(3, 2, figsize=(10, 10))  # Correct order and figure size

ax[0,0].plot(a['result'])
ax[0,0].set_title('Daylight on equator')
ax[0,0].set_ylim(0,24.5)
ax[0,0].set_xlabel('Day of Year')
ax[0,0].set_ylabel('Daylight [h]')

ax[0,1].plot(b['result'])
ax[0,1].set_title('Daylight on north pole')
ax[0,1].set_ylim(0,24.5)
ax[0,1].set_xlabel('Day of Year')
ax[0,1].set_ylabel('Daylight [h]')

ax[1,0].plot(c['result'])
ax[1,0].set_title('Daylight on south pole')
ax[1,0].set_ylim(0,24.5)
ax[1,0].set_xlabel('Day of Year')
ax[1,0].set_ylabel('Daylight [h]')

ax[1,1].plot(d['result'])
ax[1,1].set_title('Above the Arctic Circle')
ax[1,1].set_ylim(0,24.5)
ax[1,1].set_xlabel('Day of Year')
ax[1,1].set_ylabel('Daylight [h]')

ax[2,0].plot(e['result'])
ax[2,0].set_title('Close to equator')
ax[2,0].set_ylim(0,24.5)
ax[2,0].set_xlabel('Day of Year')
ax[2,0].set_ylabel('Daylight [h]')

# Remove the empty subplot (optional)
fig.delaxes(ax[2,1])  # This removes the unused subplot
plt.suptitle('Sanity Checks of the Model')
plt.xlabel('Day')
plt.ylabel('Daylight [h]')
plt.tight_layout()  # Adjust layout to prevent overlapping titles
plt.show()



script_dir = os.path.dirname(os.path.abspath(__file__))

'''
    Warsaw  - compare to date from Internet
'''
df_warsaw = pd.read_csv(script_dir + '\Warsaw_daylight.csv', header=None)
df_warsaw.columns = ['day', 'daylight']

w = plot_daylight(52.2297*np.pi/180, df=df_warsaw, num_years=1, title='Daylight in Warsaw', vis=True)
print(f"WARSAW:\nRMSE in minutes: {w['rmse']*60} \nMax difference in minutes: {w['max_diff']*60}\n")


'''
    Sydney - compare to date from Internet
'''
df_sydney = pd.read_csv(script_dir +'\Sydney_daylight.csv', header=None)
df_sydney.columns = ['day', 'daylight']
s = plot_daylight(-33.8688*np.pi/180, df=df_sydney, num_years=1, title='Daylight in Sydney', vis=True)
print(f"SYDNEY:\nRMSE in minutes: {s['rmse']*60} \nMax difference in minutes: {s['max_diff']*60}")



