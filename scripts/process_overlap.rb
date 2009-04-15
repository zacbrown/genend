#!/usr/bin/env ruby

require 'rubygems'
require 'gruff'
#require 'gsl'

SPACINGS = 50

class Record
  attr_accessor :score, :match

  def initialize(score, match)
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

  def sum
    inject( nil ) { |sum,x| sum ? sum+x : x }
  end

  def mean
    sum / size
  end
end

def process_groups(score_array, match_array)
  small = score_array.mean
  len = score_array.size
  ret_array = [[],[],[],[],[],[]]
  match_index = 0

  match_array.each { |row|
    ind = 0
    row.each do |val|
      if val == nil
        ret_array.each { |array| array << small }
        ind += 1
      else
        ret_array[ind] << score_array[match_index]
        ind += 1
      end
    end
    match_index += 1
  }

  return ret_array
end

# def poly_fit(data, degree)
#   y = Vector(data)
#   last_val = data.length - 1
#   x = Vector((0..last_val).to_a)

#   coef, err, chisq, status = MultiFit.polyfit(x, y, degree)
#   x2 = Vector.linspace(1, last_val, SPACINGS)

#   return coef.eval(x2)
# end

file = ARGV[0]
split_num = ARGV[1]
piece_size = file.split('-')[1]
date = file.split('-')[2].chomp('.dat')
header = file.split('-')[0]
NUM_KMERS = 31050

score_vals = {'3'=>[],'4'=>[],'5'=>[],'6'=>[],'7'=>[],'8'=>[]}

File.open(file, 'r').each { |line|
  tokens = line.split(' ')
  score_vals[tokens[0]] << Record.new(tokens[7].to_f, tokens[1..6])
}

graph_hash = {}
short_len = 1000000000

score_vals.each { |key,value|
  next if value == []
#  value.sort! { |a,b| a.score <=> b.score }
  puts "processing: #{key}-mers"

  score_arr = Array.new
  match_arr = Array.new
  value.each { |record|
    score_arr << record.score
    match_arr << record.match
  }

  short_len = score_arr.size if score_arr.size < short_len
  graph_hash[key] = process_groups(score_arr, match_arr)
}

colors = {:my_green =>'#00FF00',:my_pale_green => '#66CC00',
  :my_pale_yellow => '#CCFF00', :my_yellow => '#FFFF00',
  :my_orange => '#FFA500', :my_red => '#FF0000'}

# setup data
graph_hash.each { |key,value|
  graph = Gruff::Area.new('1200x800')
  graph.theme_keynote #2E37FE
  graph.replace_colors ["#ffffff", "#00ff00", "#3333FF",
                      "#ff00ff", "#ffff00", "#ff0000"]

  graph.title = "#{header} #{key}-mer 10000 bp Overlap Taxonomic Classification"
  graph.title_font_size = 20
  graph.marker_font_size = 16
  graph.x_axis_label = "relative piece location"
  graph.y_axis_label = "raw score"
  graph.maximum_value = 1.85
  graph.minimum_value = 1.55

  array_dat = value.reverse
  array_dat.each { |tmp_arr|
    tmp_arr.map! { |item| ("%5.6f" % ((item /= 10000) + 3)).to_f }
    tmp_mean = tmp_arr.mean
    tmp_arr.map! { |item|
       if item < graph.maximum_value and item > graph.minimum_value
         item
       else
         item = tmp_mean
       end
    }
  }

  graph.data "phylum", array_dat[0]
  graph.data "class", array_dat[1]
  graph.data "order", array_dat[2]
  graph.data "family", array_dat[3]
  graph.data "genus", array_dat[4]
  graph.data "species", array_dat[5]


  graph.write("result-#{header}-#{piece_size}-#{key}-#{date}.png")

}


