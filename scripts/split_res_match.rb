#!/usr/bin/env ruby

file_toks = ARGV[0].split('-')
size = file_toks[1]
file_h = nil
cur_spec = ""

File.open(ARGV[0], 'r').each { |line|
  tokens = line.split(' ')
  if cur_spec != tokens[1]
    cur_spec = tokens[1]
    file_h.close
    file_h = File.open('#{tokens[1]}-#{size}.dat', 'w')
  end

  file_h.write(line)
}

file_h.close
