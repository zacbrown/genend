#!/usr/bin/env ruby
spec_1_fh = File.open('NC_009613-10000-match.dat', 'a')
spec_2_fh = File.open('NC_004603-10000-match.dat', 'a')
spec_3_fh = File.open('NC_007622-10000-match.dat', 'a')
spec_4_fh = File.open('NC_003295-10000-match.dat', 'a')
spec_5_fh = File.open('NC_005823-10000-match.dat', 'a')

File.open(ARGV[0], 'r').each { |line|
  tokens = line.split(' ')

  puts tokens[1]
  spec_1_fh.write(line) if tokens[1] == 'NC_009613'
  spec_2_fh.write(line) if tokens[1] == 'NC_004603.2'
  spec_3_fh.write(line) if tokens[1] == 'NC_007622'
  spec_4_fh.write(line) if tokens[1] == 'NC_003295'
  spec_5_fh.write(line) if tokens[1] == 'NC_005823.2'
}

spec_1_fh.close
spec_2_fh.close
spec_3_fh.close
spec_4_fh.close
spec_5_fh.close
