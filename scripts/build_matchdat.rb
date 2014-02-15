#!/usr/bin/env ruby

file = ARGV[0]
output_h = File.open("id-"+file, 'w')

File.open(file).each { |line|
  tokens = line.split(' ')
  if tokens[1] == tokens[2]
    output_h.write(line)
  end
}

output_h.close
