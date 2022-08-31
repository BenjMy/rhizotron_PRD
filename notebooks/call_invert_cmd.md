# ---------------------- All cycle ERT analysis ----------------------
## No PRD intro protocol
python invert.py -cycle 0 1 2 3 4 5 6 7 8 9 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5

## filter_seq
python invert.py -cycle 0 1 2 3 4 5 6 7 8 9 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5
python invert.py -cycle 4 5 6 7 8 9 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5

## filter_seq_rec
python invert.py -cycle 0 1 2 3 4 5 6 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 4 5 6 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5


# ---------------------- Time lapse ERT analysis ----------------------

# time zone intro before PRD
python invert.py -cycle 0 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

# time zone 1
python invert.py -cycle -99 -startD 21/6/2022,13:50 -endD 26/6/2022,14:50 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

## before/after irr cycle5/6
python invert.py -cycle -99 -startD 21/6/2022,13:50 -endD 23/6/2022,10:00 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

## before/after irr cycle6/7
python invert.py -cycle -99 -startD  29/6/2022,9:00 -endD 30/6/2022,9:00 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5



# time zone 2
python invert.py -cycle -99 -startD 5/7/2022,13:50 -endD 9/7/2022,14:50 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5

## before/after irr cycle5/6
python invert.py -cycle -99 -startD 5/7/2022,13:50 -endD 6/7/2022,14:50 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5



# per cycle
python invert.py -cycle 0 1 2 3 4 5 6 7 8 9 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 4 5 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 6 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5



# ---------------------- ICSD analysis ----------------------
python invert.py -cycle 3 4 5 6 7 8 9 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5
python invert.py -cycle 0 1 2 3 4 5 6 7 8 9 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5
python invert.py -cycle 5 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5




# ---------------------- 1 cycle analysis ----------------------
python invert.py -cycle 8 -TL 1 -icsd 0 -reprocessed 0 -filter_seq 0 -filter_seq_rec 1 -recErr 5
python invert.py -cycle 8 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5



# ---------------------- Synthetic case ----------------------
#'%d/%m/%Y,%H:%M'

python invert_synth.py -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -scenario A
python invert_synth.py -cycle 7 -startD 29/6/2022,13:50 -endD 30/6/2022,14:50 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -scenario A
python invert_synth.py -cycle 7 -startD 29/6/2022,13:50 -endD 30/6/2022,14:50 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5 -scenario B




