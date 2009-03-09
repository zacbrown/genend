#!/usr/bin/env ruby

require 'rubygems'
require 'gruff'

class Record
  attr_accessor :score, :match

  def initialize(match, score)
    @score = score
    @match = match
  end
end

class Array
  def count_true
    ret_val = 0
    self.each { |val|
      ret_val += 1 if val == true
    }
    return ret_val
  end
end

def process_groups(score_array, match_array)
  len = score_array.size
  chunk_size = len / 150
  ret_array = Array.new

  cur_chunk_ind = 0

  while cur_chunk_ind < len
    tmp_score_array = score_array[cur_chunk_ind, chunk_size - 1]
    tmp_match_array = match_array[cur_chunk_ind, chunk_size - 1]

    count = tmp_match_array.count_true
    perc = ("%5.4f" % (count.to_f / chunk_size.to_f * 100)).to_f
    for i in cur_chunk_ind .. (cur_chunk_ind + chunk_size - 1)
      ret_array[i] = perc
    end

    cur_chunk_ind += chunk_size
  end

  return ret_array
end


file = ARGV[0]
split_num = ARGV[1]
piece_size = file.split('-')[1]
date = file.split('-')[2].chomp('.dat')
header = file.split('-')[0]

score_vals = {'3'=>[],'4'=>[],'5'=>[],'6'=>[],'7'=>[],'8'=>[]}

File.open(file, 'r').each { |line|
  tokens = line.split(' ')
  score_vals[tokens[0]] << Record.new(tokens[1] == tokens[2], tokens[3].to_f)
}

graph_hash = {}
low = 0.0; high = 0.0;

score_vals.each { |key,value|
  next if value == []
  value.sort! { |a,b| a.score <=> b.score }
  puts "processing: #{key}-mers"

  score_arr = Array.new
  match_arr = Array.new
  value.each { |record|
    score_arr << record.score
    match_arr << record.match
  }

  low = ("%5.2f" % score_arr[0]).to_f if key == '3'
  high = ("%5.2f" % score_arr[score_arr.size - 1]).to_f if key == '8'

  graph_hash[key] = process_groups(score_arr, match_arr)
}

graph = Gruff::Line.new
graph.theme_keynote
graph.title = "#{header} - #{piece_size} bp pieces - 100 draws"
graph.x_axis_label = "probability score"
graph.y_axis_label = "percent identified"

# setup data
graph_hash.each { |key,value|
  graph.data "#{key}-mer", value
}

# prep labels for x axis
incr = ("%5.2f" % ((high - low) / 5.0)).to_f
puts "low: #{low} || high: #{high}"
graph_chunks = graph_hash['6'].size / 5
graph_x_vals = (1..4).to_a
graph_x_vals.map! { |val| val * graph_chunks }

cur_score_val = low
graph.labels[0] = low.to_s
graph.labels[graph_hash['8'].size - 1] = high.to_s
graph_x_vals.each { |val|
  graph.labels[val] = (cur_score_val + incr).to_s
  cur_score_val += incr
}

graph.write("result-#{header}-#{piece_size}-#{date}.png")

