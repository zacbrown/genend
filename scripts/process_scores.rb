#!/usr/bin/env ruby

require 'rubygems'
require 'gruff'
require 'gsl'

SPACINGS = 50

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
  chunk_size = len / SPACINGS
  ret_array = Array.new

  cur_chunk_ind = 0

  while cur_chunk_ind < len
    if cur_chunk_ind + chunk_size > len
      tmp_score_array = score_array[cur_chunk_ind, len - 1]
      tmp_match_array = match_array[cur_chunk_ind, len - 1]
    else
      tmp_score_array = score_array[cur_chunk_ind, chunk_size - 1]
      tmp_match_array = match_array[cur_chunk_ind, chunk_size - 1]
    end

    count = tmp_match_array.count_true
    perc = ("%5.4f" % (count.to_f / chunk_size.to_f * 100)).to_f
    for i in cur_chunk_ind .. (cur_chunk_ind + chunk_size - 1)
      ret_array[i] = perc
    end

    cur_chunk_ind += chunk_size
  end

  return ret_array
end

def poly_fit(data, degree)
  y = Vector(data)
  last_val = data.length - 1
  x = Vector((0..last_val).to_a)

  coef, err, chisq, status = MultiFit.polyfit(x, y, degree)
  x2 = Vector.linspace(1, last_val, SPACINGS)

  return coef.eval(x2)
end

file = ARGV[0]
split_num = ARGV[1]
piece_size = file.split('-')[1]
date = file.split('-')[2].chomp('.dat')
header = file.split('-')[0]
NUM_KMERS = 31050

score_vals = {'3'=>[],'4'=>[],'5'=>[],'6'=>[],'7'=>[],'8'=>[]}

File.open(file, 'r').each { |line|
  tokens = line.split(' ')
  score_vals[tokens[0]] << Record.new(tokens[1] == tokens[2], tokens[3].to_f)
}

graph_hash = {}
low = 0.0; high = 0.0;
short_len = 1000000000

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

  short_len = score_arr.size if score_arr.size < short_len
  low = ("%5.2f" % score_arr[0]).to_f if key == '3'
  high = ("%5.2f" % score_arr[score_arr.size - 2]).to_f if key == '8'
  graph_hash[key] = process_groups(score_arr, match_arr)
}

# truncate arrays that are longer than the shortest in the group
graph_hash.each { |key,value|
  if value.size > short_len
    diff = value.size - short_len
    value.slice! 0, diff
  end
}

graph = Gruff::Line.new
graph.theme_keynote #2E37FE
graph.replace_colors ["#ff0000", "#00ff00", "#3333FF",
                      "#ff00ff", "#ffff00", "#ffffff"]
graph.title = "#{header} - #{piece_size} bp pieces - 100 draws"
graph.x_axis_label = "probability score"
graph.y_axis_label = "percent identified"
graph.maximum_value = 100
graph.minimum_value = 0

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

