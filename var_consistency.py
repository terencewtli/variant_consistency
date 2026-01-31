import argparse
import con_counts
import split_con
import extract_df

def main():
    parser = argparse.ArgumentParser(description='Pipeline designed to analyze \
                                    output of genotype-based demultiplexing methods.')
    parser.add_argument(
        'step', 
        choices=['con_counts', 'split_con', 'extract_df'],
        help='Specify which step'
    )
    args = parser.parse_args()

    if args.step == 'con_counts':
        step1.run()
    elif args.step == 'split_con':
        step2.run()
    elif args.step == 'extract_df':
        step3.run()

if __name__ == "__main__":
    main()
