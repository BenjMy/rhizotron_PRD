Anisotropy analysis
======================


python invert.py -cycle 6 7 8 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5


All cycles ERT analysis
======================


No PRD intro protocol
----------------------
python invert.py -cycle 0 1 2 3 4 5 6 7 8 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5

Filter_seq
----------------------
python invert.py -cycle 0 1 2 3 4 5 6 7 8 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5
python invert.py -cycle 4 5 6 7 8 9 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5

Filter_seq_rec
----------------------
python invert.py -cycle 0 1 2 3 4 5 6 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 4 5 6 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5



Time lapse ERT analysis
=======================

# time zone 1 --> before after irrigation end cycle 4 and to end of 5 IRRIGATION RIGHT
python invert.py -cycle -99 -startD 21/6/2022,13:50 -endD 29/6/2022,10:13 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1

python invert.py -cycle -99 -startD 21/6/2022,13:50 -endD 29/6/2022,10:13 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 2

##  time zone 2 --> before/after irr cycle5/6 IRRIGATION LEFT
python invert.py -cycle -99 -startD  29/6/2022,9:00 -endD 5/7/2022,17:00 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1

python invert.py -cycle -99 -startD  29/6/2022,9:00 -endD 5/7/2022,17:00 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1


# time zone 3 --> before after irrigation end cycle 6 and to end of 7 IRRIGATION RIGHT [6/7]
python invert.py -cycle -99 -startD 5/7/2022,16:00 -endD 11/7/2022,12:05 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle -99 -startD 5/7/2022,16:00 -endD 11/7/2022,12:05 -TL 1 -TLreg 2 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5


# time zone 4 --> before after irrigation end cycle 7 and to end of 8 IRRIGATION LEFT
python invert.py -cycle -99 -startD 11/7/2022,11:19 -endD 12/7/2022,13:05 -TL 1 -TLreg 1 -icsd 0 -reprocessed 1 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle -99 -startD 11/7/2022,11:19 -endD 12/7/2022,13:05 -TL 1 -TLreg 2 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5



# time zone intro before PRD
python invert.py -cycle -1 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1
python invert.py -cycle 0 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1
python invert.py -cycle 1 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1
python invert.py -cycle 2 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1
python invert.py -cycle 3 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5 -TLreg 1

# partial irrigation: 19/05/2022 17:00-17:30 200 ml through the first 4 upper holes (left side)
python invert.py -cycle -99 -startD 25/5/2022,13:29 -endD 1/6/2022,12:51 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5



# per cycle
python invert.py -cycle -1 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 0 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5


# cycle -1 to 0
python invert.py -cycle -99 -startD 19/5/2022,15:37 -endD 25/5/2022,13:31 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

# cycle 0 to 1
python invert.py -cycle -99 -startD 25/5/2022,13:29 -endD 1/6/2022,12:51 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

# cycle 1 to 2
python invert.py -cycle -99 -startD 1/6/2022,12:49 -endD 8/6/2022,10:01 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

# cycle 2 to 3
python invert.py -cycle -99 -startD 8/6/2022,09:59 -endD  8/6/2022,12:31 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

# cycle 3 to 4
python invert.py -cycle -99 -startD 15/6/2022,16:19 -endD  22/6/2022,16:11 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5



python invert.py -cycle -99 -startD 25/5/2022,13:29 -endD 11/7/2022,12:05 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 1 2 -TL 1 -TLreg 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 0 1 2 3 4 5 6 7 8 9 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 4 5 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 6 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5




ICSD analysis
==============

python invert.py -cycle 3 4 5 6 7 8 9 -TL 0 -icsd 1 -reprocessed 1 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -pareto 0 -wr 1
python invert.py -cycle 3 4 5 6 7 8 9 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -pareto 0 -wr 1
 
 
# ---------------------- 1 cycle analysis ----------------------
python invert.py -cycle 8 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 8 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5




Synthetic case
==============
#'%d/%m/%Y,%H:%M'

python invert_synth.py -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -scenario A

python invert_synth.py -cycle 7 -startD 29/6/2022,13:50 -endD 30/6/2022,14:50 -reprocessed 1 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -scenario B


Plot papers ERT PRD effects + Archie
============================
python invert.py -cycle 6 7 8 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -petro 1


TEST
====

python invert.py -cycle 7 -TL 0 -icsd 1 -reprocessed 1 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -startD 29/6/2022,14:14 -endD 29/6/2022,15:03 -pareto 0 -wr 10
 




