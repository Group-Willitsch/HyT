# -*- coding: utf-8 -*-
"""
Classfile for calculation of the trigger sequence from the delay tab in the program

21.03.2019: Tested on experiment, working but not all cases are solved.
            changing in standard triggering scheme 3b -> 10b causes an error
            the order of solving the puzzle has to be determined before starting the solving

            TBD

22.03.2019: A background option is implemented. One gives a channel number for the BG, the sequence is solved
            and in the end a channel on (1) is inverted (0) or vice versa.

            TBD

            In the decelerator sequence it should check for other and give a minimal distance of a few us to prevent the error
            instruction delay too small from the pulseblaster

            TBD


"""
import numpy as np
import readfile as rf
import pandas as pd
import pulseblaster as pb
from treelib import Node, Tree

units = {'ns': 1e-9, 'us': 1e-6, 'ms': 1e-3, 's': 1}


def indexlist(lists, element):
    '''
    creates a list of indices with a certain value from given list
    '''
    indices = [i for i, x in enumerate(lists) if x == element]
    return indices


class sequence():
    def __init__(self, units1, chno1, time1, units2, chno2, time2, freqs,
                 pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, bgchannel,
                 scantime, t2jumpfile):

        self.units1 = units1
        self.chno1 = chno1
        self.time1 = time1
        self.units2 = units2
        self.chno2 = chno2
        self.time2 = time2
        self.freqs = freqs
        self.pulsenature = pulsenature
        self.checkboxvals = checkboxvals
        self.posnegvals = posnegvals
        self.scanchannel = scanchannel
        self.relchannel = relchannel
        self.bgchannel = bgchannel
        self.scantime = scantime
        self.t2jumpfile = t2jumpfile

        self.freqs = freqs
        self.freqmin = min(freqs)
        self.freqmax = max(freqs)
        self.pulsetime = 1 / self.freqmin  # the pulse blaster has to be programmed for this amount of time to get the reprate
        # This means that from the last pulse, the sequence has to be filled with nothing for some time

        self.events = {}  # create a dictionary for the events such as '100b': 400e-6 which would mean channel 1, event 00, begining at 400 us.

        self.ch_done = []  # list with channels which were already treated
        self.lastindex = len(self.chno1)
        self.channels_in_use = [index for index in range(self.lastindex) if self.checkboxvals[index]]
        self.ch_todo = [index for index in range(self.lastindex) if self.checkboxvals[index]]

        self.starttag = 'T0'  # begins by searching a starttag T0 in chno1
        self.t0_index = None  # index of T0 in chno1 (first channel number in trigger table)

    def find_T0(self):
        # finds the tag T0 as starting point
        # returns True if succesful, False otherwise
        indexes = indexlist(self.chno1, self.starttag)
        checked = [index for index in indexes if self.checkboxvals[index]]

        if len(checked) > 1:
            print('several T0 checked')
            return False
        elif not checked:
            print('No T0 checked')
            return False

        self.t0_index = checked[0]
        return True

    def check_freqs(self):
        # Check if the times set are smaller than the repetition time (1/freq)
        badchannels = []  # channels with problem

        for i in self.channels_in_use:
            if (self.time1[i] >= 1 / self.freqs[i]) or (self.time2[i] >= 1 / self.freqs[i]):
                badchannels.append(i + 1)
        return badchannels

        # for i in range(len(self.time1)):
        #     if (self.time1[i] > 1/self.freqs[i]) or (self.time2[i] > 1/self.freqs[i]):
        #         badchannels.append(i+1)

        # if badchannels:
        # print ('reprate: bad channels are %s'%badchannels)
        # return False, badchannels
        # else:
        # return True, None

    def t0_def(self):
        runs = int(self.freqs[self.t0_index] / self.freqmin)

        delay_b = self.time1[self.t0_index]
        delay_e = self.time2[self.t0_index]
        n = 0

        channel = self.t0_index + 1

        if self.scanchannel in '%ib' % (self.t0_index + 1):
            delay_b = self.scantime
            channel = self.relchannel[:-1]

        elif self.scanchannel in '%ie' % (self.t0_index + 1):
            delay_e = self.scantime
            channel = self.relchannel[:-1]

        while n < runs:
            if self.pulsenature[self.t0_index] == 'normal pulse':
                n += 1
                self.events['%ib%03i' % (channel, n)] = (n - 1) * 1 / self.freqs[self.t0_index] + delay_b
                # event channel begin with nr %03i (3 digit integer) e.g. ch1b001, ch1b002 etc,...

                self.events['%ie%03i' % (channel, n)] = self.events['%ib%03i' % (self.t0_index + 1, n)] + delay_e
                # end of the channel e.g. ch1e001 etc.
            n += 1

        self.ch_done.append(self.t0_index)
        self.ch_todo.remove(self.t0_index)
        return

    def tag_and_time(self, index):
        laserchan = '3b'
        thresh = 15e-6
        # defines the begining and end of a channel based on the trigger table
        if not self.checkboxvals[index]:  # basically not necessary, since ch_todo should only include checked channels
            self.ch_done.append(index)
            return False

        delay_b = round(self.time1[index], 7)
        delay_e = round(self.time2[index], 7)

        ch_b = self.chno1[index]
        ch_e = self.chno2[index]

        # XXX

        #        if self.scanchannel in '%ib'%(index+1):
        ##            print('scan b')
        #            delay_b = self.scantime
        #            ch_b = self.relchannel
        #
        #        elif self.scanchannel in '%ie'%(index+1):
        ##            print('scan e')
        #            delay_e = self.scantime
        #            ch_e = self.relchannel

        runs = int(self.freqs[index] / self.freqmin)
        # number of runs (e.g. ch1 reprate 1 Hz, rest 10 Hz -> sequence is 1 s long and ch1 is triggered once, the rest 10 times)
        n = 0  # index for runs

        if self.pulsenature[index] == 'normal pulse':
            m = 0

            while n < runs:
                if n == 0:
                    n += 1
                    self.ch_done.append(index)
                    self.ch_todo.remove(index)

                    self.events['%ib%03i' % (index + 1, n)] = self.events['%s%03i' % (ch_b, n)] + delay_b
                    self.events['%ie%03i' % (index + 1, n)] = self.events['%s%03i' % (ch_e, n)] + delay_e

                    if self.events['%ib%03i' % (index + 1, n)] >= 1 / self.freqs[index]:
                        m = np.floor(self.events['%ib%03i' % (index + 1, n)] * self.freqs[index])

                    elif self.events['%ie%03i' % (index + 1, n)] >= 1 / self.freqs[index]:
                        m = np.floor(self.events['%ie%03i' % (index + 1, n)] * self.freqs[index])

                    # on the first pulse check if freq>fremin and
                    # the first pulse is defined relative to something later than 1/freq (e.g. flashlamp 10 Hz before Qsw 2 Hz)
                    self.events['%ib%03i' % (index + 1, n)] -= m * 1 / self.freqs[index]
                    self.events['%ie%03i' % (index + 1, n)] -= m * 1 / self.freqs[index]


                else:
                    n += 1
                    n2 = n - m

                    self.events['%ib%03i' % (index + 1, n)] = self.events['%ib%03i' % (index + 1, n - 1)] + 1 / \
                                                              self.freqs[index]
                    self.events['%ie%03i' % (index + 1, n)] = self.events['%ie%03i' % (index + 1, n - 1)] + 1 / \
                                                              self.freqs[index]

        elif 'decelerator' in self.pulsenature[index]:
            n_dec = 0
            m = 0

            ####Select if you want to read the fortran generated sequence (classic) or python generated sequence
            classic = True
            if 'H+' in self.pulsenature[index]:
                if classic:
                    dectime, decdur = rf.readdecfile(self.t2jumpfile, 'H+')
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/v_time_470',dectime)
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/v_duration_470',decdur)
                else:
                    print("Note Python sequence, be carefule")
                    dectime = np.load('./input/input_test_py_code/v_time_425.npy', allow_pickle=True)
                    decdur = np.load('./input/input_test_py_code/v_duration_425.npy')

            if 'H-' in self.pulsenature[index]:
                if classic:
                    dectime, decdur = rf.readdecfile(self.t2jumpfile, 'H-')
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/h_time_470',dectime)
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/h_duration_470',decdur)
                else:
                    print("Note Python sequence, be carefule")
                    dectime = np.load('./input/input_test_py_code/h_time_425.npy', allow_pickle=True)
                    decdur = np.load('./input/input_test_py_code/h_duration_425.npy')
            if 'V+' in self.pulsenature[index]:
                if classic:
                    dectime, decdur = rf.readdecfile(self.t2jumpfile, 'V+')
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/v_time_470',dectime)
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/v_duration_470',decdur)
                else:
                    print("Note Python sequence, be carefule")
                    dectime = np.load('./input/input_test_py_code/v_time_425.npy', allow_pickle=True)
                    decdur = np.load('./input/input_test_py_code/v_duration_425.npy')

            if 'V-' in self.pulsenature[index]:
                if classic:
                    dectime, decdur = rf.readdecfile(self.t2jumpfile, 'V-')
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/h_time_470',dectime)
                    # np.save('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/0deg_470/h_duration_470',decdur)
                else:
                    print("Note Python sequence, be carefule")
                    dectime = np.load('./input/input_test_py_code/h_time_425.npy', allow_pickle=True)
                    decdur = np.load('./input/input_test_py_code/h_duration_425.npy')

            while n < runs:
                if n == 0:

                    n += 1
                    self.ch_done.append(index)
                    self.ch_todo.remove(index)

                    self.events['%ib%03i' % (index + 1, n)] = self.events['%s%03i' % (ch_b, n)] + delay_b

                    if self.events['%ib%03i' % (index + 1, n)] >= 1 / self.freqs[index]:
                        m = np.floor(self.events['%ib%03i' % (index + 1, n)] * self.freqs[index])

                    n2 = n - m
                    n_dec += 1

                    for i in range(dectime.shape[0]):
                        if self.events['%s%03i' % (self.chno1[index], 1)] + self.time1[index] + (n2 - 1) * 1 / \
                                self.freqs[index] + dectime[i] < self.events['3b001'] - thresh:
                            self.events['%ib%03i' % (index + 1, n_dec)] = self.events[
                                                                              '%s%03i' % (self.chno1[index], 1)] + \
                                                                          self.time1[index] + (n2 - 1) * 1 / self.freqs[
                                                                              index] + dectime[i]
                            self.events['%ie%03i' % (index + 1, n_dec)] = self.events['%ib%03i' % (index + 1, n_dec)] + \
                                                                          decdur[i]
                            n_dec += 1

                else:
                    n2 = n - m
                    n_dec += 1

                    for i in range(dectime.shape[0]):
                        #                        if self.events['%s%03i'%(self.chno1[index],1)]+self.time1[index]+(n2-1)*1/self.freqs[index]+dectime[i]
                        if self.events['%s%03i' % (self.chno1[index], 1)] + self.time1[index] + (n2 - 1) * 1 / \
                                self.freqs[index] + dectime[i] < self.events['3b001'] - thresh:
                            self.events['%ib%03i' % (index + 1, n_dec)] = self.events[
                                                                              '%s%03i' % (self.chno1[index], 1)] + \
                                                                          self.time1[index] + (n2 - 1) * 1 / self.freqs[
                                                                              index] + dectime[i]
                            self.events['%ie%03i' % (index + 1, n_dec)] = self.events['%ib%03i' % (index + 1, n_dec)] + \
                                                                          decdur[i]
                            n_dec += 1
        #                    print(self.pulsenature[index], n_dec)
        return True

    def format_pulseblaster(self):
        # we have to transform the events into the format, that the pulse blaster accepts e.g. channels 1 and 3 on: 1010000 for a duration.
        # the sequence also has to be finished with 000000 or 010000 if 1 is a negative pulse
        data = pd.DataFrame([self.events.values(), self.events.keys()]).T.sort_values(by=0).reset_index()
        #        print(data)

        digits = 7  # 1 ns resolution
        data[0] = np.around(np.array(data[0], dtype=np.float64), decimals=digits)
        data = data.groupby([0])[1].apply(list)

        #        print(data)

        #    Correct until here
        # ~ create the pulses for all channels
        '''
            extract int from string:
            [int(s) for s in str.split() if s.isdigit()]
        '''

        activechannels = []  # list for the channels active during one event
        self.output = []  # list for output: channelflag and duration

        times = np.array(data.keys())
        times = times - np.min(times)  # time starts at zero
        lastdig = 4

        #        print(times)

        number = 0
        for events in data:
            if events == list(data)[-1]:
                continue

            duration = times[number + 1] - times[number]
            number += 1

            for keys in events:
                channel = int(keys[:-lastdig])
                # print(keys)
                flag = keys[-lastdig]
                m = keys[-lastdig + 1:]

                if flag == 'b':
                    activechannels.append(channel)

                if flag == 'e':
                    activechannels.remove(channel)

            chflag = pb.channelflag(activechannels)
            for chno in range(len(self.posnegvals)):
                if self.checkboxvals[chno]:
                    ch = list(chflag)

                    if self.posnegvals[chno] == 'Negative':
                        if chflag[chno] == '1':
                            ch[chno] = '0'
                        else:
                            ch[chno] = '1'

                        if chno + 1 == self.bgchannel:  # check for background channel, invert
                            #                            print('here 1')
                            ch[chno] = '1'

                    else:  # if not negative
                        if chno + 1 == self.bgchannel:  # check for background channel, invert
                            #                            print('here bg', chno)

                            ch[chno] = '0'
                    chflag = "".join(ch)

            if self.checkboxvals[channel - 1]:
                self.output.append([chflag, float(duration)])

        '''
    #        finish whole sequence and fill with waittime
        '''
        chflag = pb.channelflag([])
        for chno in range(len(self.posnegvals)):
            if self.checkboxvals[chno]:
                ch = list(chflag)

                if self.posnegvals[chno] == 'Negative':
                    if chflag[chno] == '1':
                        ch[chno] = '0'

                    else:
                        ch[chno] = '1'

                    if chno + 1 == self.bgchannel:  # check for background channel, invert
                        ch[chno] = '1'
                #                            print('here 3')

                else:  # if not negative
                    if chno + 1 == self.bgchannel:  # check for background channel, invert
                        ch[chno] = '0'
                #                            print('here 4')

                chflag = "".join(ch)

        duration = 1 / self.freqmin - np.max(times)
        self.output.append([chflag, float(duration)])
        self.output = np.array(self.output).T
        #        print(self.output)

        return

    def seq(self):
        # calculates the sequence using the input form the gui a scan channel can be specified relative to another channel (e.g. 1e rel 1b)
        if not self.find_T0():  # check if succesful in finding T0
            return False
        badchannels = self.check_freqs()  # check if times set are smaller than the repetition time
        if badchannels:
            raise Exception("Error: Trigger frequency too high!")

        self.t0_def()
        index = self.t0_index  # start with channels in relation to t0 channel
        # print('here'+str(index))

        self.find_channel_relation()
        for key, value in self.channel_relation_dict.items():
            for v in value:
                if v in self.channel_relation_dict and v != 5:
                    try:
                        self.tag_and_time(v)
                    except:
                        print("tag_and_time(%i) error" % v)
                    # self.channel_relation_tree.create_node(identifier=v, parent= key)


        # while len(self.ch_done) < len(self.channels_in_use):
        #     try:
        #         index = indexlist(self.chno1, '%ib' % (index + 1))[
        #             0]  # look in chno1 list for a channel that has his beginning relative
        #         # to the defined channel (first time will be relative to the channel with T0 (channel numbers are index+1 !!!)
        #         # ~ take first in this list
        #     except IndexError:
        #         print('Index Error channel flags: %i' % (index + 1))
        #         # the loop breaks
        #         break
        #
        #     try:
        #         print(index)
        #         self.tag_and_time(index)
        #     except:
        #         print("tag_and_time(%i) error" % index)
        todo = [item for item in self.ch_todo]  # channels which have to be calculated
        print(todo)

        # the rest of the channels which were not calculated so far are done now
        for index in todo:
            self.tag_and_time(index)

        # we have to transform the events now into the format, that the pulse blaster accepts e.g. channels 1 and 3 on: 1010000 for a duration.
        # the sequence also has to be finished with 000000 or 010000 if 1 is a negative pulse
        #        print(self.events)
        try:
            self.format_pulseblaster()
        except:
            print("format_pulseblaster() error")

        return True

    def find_channel_relation(self):
        self.channel_relation_dict = {}  # key: channel_index, e.g. 1b, value: a list of channel indices that depends on the channel specified by the key.
        for ch in self.channels_in_use:
            chs = [i for i, x in enumerate(self.chno1) if str(ch + 1) == x[:-1]]
            if chs:
                self.channel_relation_dict[ch] = chs
        # return self.channel_relation_dict
        print(self.channel_relation_dict)

    def create_channel_relation_tree(self):

        self.channel_relation_tree = Tree()
        self.channel_relation_tree.create_node(identifier=0, data=0)
        # self.channel_relation_tree.create_node(data=3, parent=0, identifier=3)
        # self.channel_relation_tree.show()

        for key, value in self.channel_relation_dict.items():
            for v in value:
                if v in self.channel_relation_dict:
                    self.channel_relation_tree.create_node(identifier=v, parent= key)
        self.channel_relation_tree.show()

        # num_trees = 1
        # for key, value in self.channel_relation_dict.items():
        #     exec("ch_tree_%i = Tree()" % key)
        #     exec("ch_tree_%i.create_node(data = %i, identifier = %i)" % (key, key, key))
        #     num_trees += 1
        #     for ch in value:
        #         exec("ch_tree_%i.create_node(data = %i , parent = %i, identifier = %i)" % (key, ch, key, ch))
        #     exec("ch_tree_%i.show()" % key)
        #
        # while num_trees > 1:
        #     print(num_trees)
        #     for key in self.channel_relation_dict:
        #
        #         num_trees = num_trees - 1
        #         # if key != 0 and exec("ch_tree_%i.contains(%i)" % (0, key)):
        #         #     exec("ch_tree_0.paste(%i, ch_tree_%i)" % (key, key))
        #
        # exec("ch_tree_0.show()")


if __name__ == "__main__":

    filename = './input/T2jump.dat'  # inputfile
    t2jumpfile = './input/sequencesFM_vel_12p5kV_124switches/dec_450_92/T2jump.out'
    # t2jumpfile = './input/sequences_new/deceleration_new/dec_450_62p5_30/T2jump.out' # decelerator sequence file from fortran trajectory simulation

    units1 = ['us', 'us', 'us', 'us', 'us', 'us', 'us', 'us', 'us', 'us']
    units2 = units1
    chno1 = ['T0', '3b', '10b', '1b', '3b', '4b', '6b', '6b', '6b', '4b']
    chno2 = ['1b', '2b', '3b', '4b', '5b', '3b', '3b', '3b', '3b', '10b']
    time1 = [400e-6, -170e-6, 1800e-6, 130e-6, -5e-6, 500e-6, 0.0, 0.0, 0.0, 2000e-6]
    time2 = [100e-6, 10e-6, 10e-6, 70e-6, 10e-6, -5e-6, -5e-6, -5e-6, -5e-6, 100e-6]
    freqs = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
    posnegvals = ['Positive', 'Positive', 'Positive', 'Positive', 'Positive', 'Positive', 'Positive', 'Positive',
                  'Positive', 'Positive']

    scanchannel = 'No'  # '3b'#'2e'
    relchannel = 'No'  # '4b'
    bgchannel = 'No'  # 2

    checkboxvals = [True, True, True, True, True, True, True, True, True, True]
    pulsenature = ['normal pulse', 'normal pulse', 'normal pulse', 'normal pulse', 'normal pulse',
                   'burst unit decelerator (H+)', 'burst unit decelerator (H-)', 'burst unit decelerator (V+)',
                   'burst unit decelerator (V-)', 'normal pulse']

    scantime = 0

    seq = sequence(units1, chno1, time1, units2, chno2, time2, freqs,
                   pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, bgchannel,
                   scantime, t2jumpfile)

    seq.seq()
    sequencedata = seq.output
    seq.find_channel_relation()
    seq.create_channel_relation_tree()
    #    print(sequencedata)

    import plotsequence as plotseq
    from matplotlib import pyplot as plt

    xvals, yvals, newlabels = plotseq.plotsequence(sequencedata, checkboxvals)

    fig, ax = plt.subplots()
    ticks = [0.25 + i for i in range(len(xvals))]
    ax.set_yticks(ticks)
    ax.set_yticklabels(newlabels)
    plt.ylabel('Channel no.')
    plt.xlabel('time [$\mu$s]')

    no = 0
    for i in range(len(xvals)):
        plt.plot(np.array(xvals[i]) * 1e6, (np.array(yvals[i]) / 2) - (newlabels[i] - 1) + i)
    plt.show()
